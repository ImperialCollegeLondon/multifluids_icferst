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
    use solvers
    use sparse_tools
    use sparsity_patterns
    use global_parameters, only: timestep

    use field_derivatives, only : strain_rate
    
    implicit none
    
    type rheology_type
       integer :: rheology_model
       type(power_law_rheology_type), pointer :: power_law_ptr=>null()
       type(carreau_rheology_type), pointer :: carreau_ptr=>null()
       type(viscoelastic_rheology_type), pointer :: viscoelastic_ptr=>null()
    end type rheology_type

    type power_law_rheology_type
       real :: yield_stress
       real :: exponent
       real :: consistency_index
       real :: tolerance
       real :: lower_bound,upper_bound
    end type power_law_rheology_type

    type carreau_rheology_type
       real :: mu_inf
       real :: mu_0
       real :: lambda
       real :: exponent
    end type carreau_rheology_type

    type viscoelastic_rheology_type
       real :: mu_s
       real :: mu_p
       real :: lambda
    end type viscoelastic_rheology_type


    private
  
    public :: rheology_type, initialize_rheologies, finalize_rheologies, calculate_rheologies, update_old_full_stress



    contains

 subroutine project_field_full_scalar(from_field, to_field, X)
        type(scalar_field) :: from_field, to_field
        type(vector_field) :: X

        type(csr_sparsity) :: mass_sparsity
        type(csr_matrix) :: mass
        type(scalar_field) :: rhs_field
        integer :: ele
        if(from_field%mesh==to_field%mesh) then
       
       call set(to_field, from_field)
       
    else

       if (to_field%mesh%continuity<0) then
          ! DG case

          do ele=1,element_count(to_field)
             call dg_projection_ele(ele, from_field, to_field, X)
          end do

       else
          ! CG case
          call allocate(rhs_field,mesh=to_field%mesh,&
               name='ProjectionRHS')
          call project_field(from_field,to_field,X)
          call zero(rhs_field)
          mass_sparsity = make_sparsity(to_field%mesh, to_field%mesh, name="MassMatrixSparsity")
          call allocate(mass,mass_sparsity,name="MassMatrix")
          call zero(mass)
          call compute_mass(X, to_field%mesh, mass)
        
           do ele=1,element_count(to_field)
             call cg_projection_ele(ele, from_field, rhs_field, X)
          end do

          call set_solver_options(rhs_field,ksptype='cg',pctype='mg', max_its=10000,atol=1.0e-12)

          call petsc_solve(to_field,mass, rhs_field,rhs_field%option_path)

          call deallocate(rhs_field)
          call deallocate(mass)
          call deallocate(mass_sparsity)
       end if
    end if

  contains
    
    subroutine dg_projection_ele(ele, from_field, to_field, X)
      integer :: ele
      type(scalar_field), intent(in) :: from_field
      type(scalar_field), intent(inout) :: to_field    
      type(vector_field), intent(in) :: X
      
      real, dimension(ele_loc(to_field,ele), ele_loc(to_field,ele)) :: mass
      real, dimension(ele_ngi(to_field,ele)) :: detwei
      type(element_type), pointer :: to_shape

  
      
      call transform_to_physical(X, ele, detwei)
      
      to_shape=>ele_shape(to_field, ele)
      
      mass=shape_shape(to_shape, to_shape, detwei) 
      
      call invert(mass)
      
      call set(to_field,ele_nodes(to_field, ele), &
           matmul(mass, &
           shape_rhs(to_shape, ele_val_at_quad(from_field,ele)*detwei)))

    end subroutine dg_projection_ele
    
    subroutine cg_projection_ele(ele, from_field, to_field, X)
      integer :: ele
      type(scalar_field), intent(in) :: from_field
      type(scalar_field), intent(inout) :: to_field
      type(vector_field), intent(in) :: X
      
      real, dimension(ele_ngi(to_field,ele)) :: detwei
      type(element_type), pointer :: to_shape
         
      to_shape=>ele_shape(to_field, ele)
      
      call transform_to_physical(X, ele, detwei)

      call addto(to_field, ele_nodes(to_field, ele), &
           shape_rhs(to_shape, ele_val_at_quad(from_field, ele)*detwei))

    end subroutine cg_projection_ele
       
  end subroutine project_field_full_scalar

      subroutine project_field_full_vector(from_field, to_field, X)
        type(vector_field) :: from_field, to_field
        type(vector_field) :: X

        type(csr_sparsity) :: mass_sparsity
        type(csr_matrix) :: mass
        type(vector_field) :: rhs_field
        integer :: ele
        if(from_field%mesh==to_field%mesh) then
       
       call set(to_field, from_field)
       
    else

       if (to_field%mesh%continuity<0) then
          ! DG case

          do ele=1,element_count(to_field)
             call dg_projection_ele(ele, from_field, to_field, X)
          end do

       else
          ! CG case

          

          call allocate(rhs_field,mesh=to_field%mesh,&
               name='ProjectionRHS',dim=to_field%dim)
          call zero(to_field)
          call zero(rhs_field)
          mass_sparsity = make_sparsity(to_field%mesh, to_field%mesh, name="MassMatrixSparsity")
          call allocate(mass,mass_sparsity,name="MassMatrix")
          call zero(mass)
          call compute_mass(X, to_field%mesh, mass)


          if (any(from_field%val/=from_field%val)) then
             print*, 'arghh!', from_field%val
          end if
        
           do ele=1,element_count(to_field)
             call cg_projection_ele(ele, from_field, rhs_field, X)
          end do
          
          call set_solver_options(rhs_field)


          call petsc_solve(to_field,mass, rhs_field,rhs_field%option_path)

          call deallocate(rhs_field)
          call deallocate(mass)
          call deallocate(mass_sparsity)
       end if
    end if

  contains
    
    subroutine dg_projection_ele(ele, from_field, to_field, X)
      integer :: ele
      type(vector_field), intent(in) :: from_field
      type(vector_field), intent(inout) :: to_field    
      type(vector_field), intent(in) :: X
      
      real, dimension(ele_loc(to_field,ele), ele_loc(to_field,ele)) :: mass
      real, dimension(ele_ngi(to_field,ele)) :: detwei
      type(element_type), pointer :: to_shape

      integer :: dim1
      
      call transform_to_physical(X, ele, detwei)
      
      to_shape=>ele_shape(to_field, ele)
      
      mass=shape_shape(to_shape, to_shape, detwei) 
      
      call invert(mass)
      
      do dim1=1,to_field%dim
         call set(to_field, dim1,ele_nodes(to_field, ele), &
              matmul(mass, &
              shape_rhs(to_shape, ele_val_at_quad(from_field,ele,dim1)*detwei)))
      end do

    end subroutine dg_projection_ele
    
    subroutine cg_projection_ele(ele, from_field, to_field, X)
      integer :: ele
      type(vector_field), intent(in) :: from_field
      type(vector_field), intent(inout) :: to_field
      type(vector_field), intent(in) :: X
      
      real, dimension(ele_ngi(to_field,ele)) :: detwei
      type(element_type), pointer :: to_shape
         
      to_shape=>ele_shape(to_field, ele)
      
      call transform_to_physical(X, ele, detwei)

      if (any( ele_val_at_quad(from_field, ele) /= ele_val_at_quad(from_field, ele)) .or.&
           any( detwei /=detwei)) then
         print*, ele
         print*, 'detwei', detwei
         print*, 'eleval',ele_val_at_quad(from_field, ele) /= ele_val_at_quad(from_field, ele)
         print*, 'eleval',ele_val_at_quad(from_field, ele)
      end if

         

      call addto(to_field, ele_nodes(to_field, ele), &
           shape_vector_rhs(to_shape, ele_val_at_quad(from_field, ele), detwei))

    end subroutine cg_projection_ele
       
  end subroutine project_field_full_vector

      subroutine project_field_full_tensor(from_field, to_field, X)
        type(tensor_field) :: from_field, to_field
        type(vector_field) :: X

        type(csr_matrix) :: mass
        type(tensor_field) :: rhs_field
        integer :: ele
        if(from_field%mesh==to_field%mesh) then
       
       call set(to_field, from_field)
       
    else

       if (to_field%mesh%continuity<0) then
          ! DG case

          do ele=1,element_count(to_field)
             call dg_projection_ele(ele, from_field, to_field, X)
          end do

       else
          ! CG case
          call allocate(rhs_field,from_field%mesh,&
               'ProjectionRHS',dim=from_field%dim)
          call zero(to_field)
          call zero(rhs_field)
          call compute_mass(X, from_field%mesh, mass)

        
           do ele=1,element_count(to_field)
             call cg_projection_ele(ele, from_field, rhs_field, X)
          end do

          call petsc_solve(to_field,mass, rhs_field)

          call deallocate(rhs_field)
          call deallocate(mass)
       end if
    end if

  contains
    
    subroutine dg_projection_ele(ele, from_field, to_field, X)
      integer :: ele
      type(tensor_field), intent(in) :: from_field
      type(tensor_field), intent(inout) :: to_field    
      type(vector_field), intent(in) :: X
      
      real, dimension(ele_loc(to_field,ele), ele_loc(to_field,ele)) :: mass
      real, dimension(ele_ngi(to_field,ele)) :: detwei
      type(element_type), pointer :: to_shape

      integer :: dim2, dim1
      
      call transform_to_physical(X, ele, detwei)
      
      to_shape=>ele_shape(to_field, ele)
      
      mass=shape_shape(to_shape, to_shape, detwei) 
      
      call invert(mass)
      
      do dim2=1,to_field%dim(2)
         do dim1=1,to_field%dim(1)
            call set(to_field, dim1, dim2,ele_nodes(to_field, ele), &
                 matmul(mass, &
                 shape_rhs(to_shape, ele_val_at_quad(from_field,ele,dim1,dim2)*detwei)))
         end do
      end do

    end subroutine dg_projection_ele
    
    subroutine cg_projection_ele(ele, from_field, to_field, X)
      integer :: ele
      type(tensor_field), intent(in) :: from_field
      type(tensor_field), intent(inout) :: to_field
      type(vector_field), intent(in) :: X
      
      real, dimension(ele_ngi(to_field,ele)) :: detwei
      type(element_type), pointer :: to_shape
         
      to_shape=>ele_shape(to_field, ele)
      
      call transform_to_physical(X, ele, detwei)

      call addto(to_field, ele_nodes(to_field, ele), &
           shape_tensor_rhs(to_shape, ele_val_at_quad(from_field, ele), detwei))

    end subroutine cg_projection_ele
       
  end subroutine project_field_full_tensor

     subroutine nullify_rheology(rheology)
        type(rheology_type), intent(inout) :: rheology

        nullify(rheology%power_law_ptr)
        nullify(rheology%carreau_ptr)
        nullify(rheology%viscoelastic_ptr)
        
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
        else if (have_option(option_path//'/rheology/non_newtonian::Carreau')) then
           rheology%rheology_model=3
           allocate(rheology%carreau_ptr)
           call get_option(option_path//&
                '/rheology/non_newtonian::Carreau/viscosity_at_zero_shear',&
                rheology%carreau_ptr%mu_0, default=0.0)
           call get_option(option_path//&
                '/rheology/non_newtonian::Carreau/viscosity_at_infinite_shear',&
                rheology%carreau_ptr%mu_inf, default=0.0)
           call get_option(option_path//&
                '/rheology/non_newtonian::Carreau/relaxation_time',&
                rheology%carreau_ptr%lambda, default=0.0)
            call get_option(option_path//&
                '/rheology/non_newtonian::Carreau/exponent',&
                rheology%carreau_ptr%exponent, default=1.0)
         else if (have_option(option_path//'/rheology/viscoelastic::Local')) then
           rheology%rheology_model=4
           allocate(rheology%viscoelastic_ptr)
           call get_option(option_path//&
                '/rheology/viscoelastic::Local/solvent_viscosity',&
                rheology%viscoelastic_ptr%mu_s, default=0.0)
           call get_option(option_path//&
                '/rheology/viscoelastic::Local/polymer_viscosity',&
                rheology%viscoelastic_ptr%mu_p, default=0.0)
           call get_option(option_path//&
                '/rheology/non_newtonian::Carreau/relaxation_time',&
                rheology%viscoelastic_ptr%lambda, default=1.0)
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
           else if (rheology(i)%rheology_model==3) then
              call allocate_and_insert_tensor_field(trim(state(i)%option_path)&
                   //'/rheology/non_newtonian[0]/tensor_field',&
                   state(i))
           else if (rheology(i)%rheology_model==4) then
              call allocate_and_insert_tensor_field(trim(state(i)%option_path)&
                   //'/rheology/viscoelastic[0]/tensor_field',&
                   state(i))
              call allocate_and_insert_tensor_field(trim(state(i)%option_path)&
                   //'/rheology/viscoelastic[1]/tensor_field',&
                   state(i))
           end if

        end do
           
      end subroutine initialize_rheologies

      subroutine finalize_rheology(rheology)
        type(rheology_type) :: rheology

        if (associated(rheology%power_law_ptr)) &
             deallocate(rheology%power_law_ptr)
        if (associated(rheology%carreau_ptr)) &
             deallocate(rheology%carreau_ptr)
        call nullify_rheology(rheology)

      endsubroutine finalize_rheology

      subroutine finalize_rheologies(rheology)
        type(rheology_type), dimension(:) :: rheology

        integer :: i

        do i=1,size(rheology)
           call finalize_rheology(rheology(i))
        end do
      end subroutine finalize_rheologies


      subroutine calculate_power_law_rheology(positions,state,vstate,viscosity,rheology,&
           continuous_velocity_mesh,material_mesh)
        type(tensor_field)  :: viscosity
        type(power_law_rheology_type) :: rheology
        type(state_type) :: state, vstate
        type(vector_field) :: positions
        type(mesh_type) :: continuous_velocity_mesh, material_mesh

        type(vector_field), pointer :: velocity
        type(tensor_field) :: strain_rate_tensor
        type(vector_field) :: continuous_velocity
        type(scalar_field) :: strain_rate_invariant
        type(scalar_field) :: scalar_viscosity, projected_scalar_viscosity

        integer :: i, j, dim


        velocity=>extract_vector_field(vstate,"Velocity")

        dim=mesh_dim(velocity)
        call allocate(continuous_velocity,mesh=continuous_velocity_mesh,&
             name="ContinuousVelocity",dim=dim)
        call allocate(strain_rate_tensor,material_mesh,name="StrainRateTensor",&
             dim=[dim,dim])
        call allocate(strain_rate_invariant,material_mesh,"StrainRateInvariant")
        call allocate(scalar_viscosity,material_mesh,"ScalarViscosity")
        call allocate(projected_scalar_viscosity,viscosity%mesh,&
             "ProjectedScalarViscosity")
        
        call project_field_full_Vector(velocity, continuous_velocity, positions)
        call strain_rate(continuous_velocity,positions,strain_rate_tensor)
        call tensor_second_invariant(strain_rate_tensor, strain_rate_invariant)
        
        call bound(strain_rate_invariant,rheology%tolerance,huge(rheology%tolerance))
        scalar_viscosity%val=(2.0*sqrt(strain_rate_invariant%val))**&
             (rheology%exponent-1.0d0)

        scalar_viscosity%val=&
             rheology%consistency_index*scalar_viscosity%val&
             +rheology%yield_stress/(2.0*sqrt(strain_rate_invariant%val))
        
        call bound(scalar_viscosity,&
             rheology%lower_bound,&
             rheology%upper_bound)


        call project_field_full_scalar(scalar_viscosity,&
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


      subroutine calculate_carreau_rheology(positions,state,viscosity,rheology,&
           continuous_velocity_mesh,material_mesh)
        type(tensor_field)  :: viscosity
        type(carreau_rheology_type) :: rheology
        type(state_type) :: state
        type(vector_field) :: positions
        type(mesh_type) :: continuous_velocity_mesh, material_mesh

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
        call allocate(strain_rate_tensor,material_mesh,name="StrainRateTensor",&
             dim=[dim,dim])
        call allocate(strain_rate_invariant,material_mesh,"StrainRateInvariant")
        call allocate(scalar_viscosity,material_mesh,"ScalarViscosity")
        call allocate(projected_scalar_viscosity,viscosity%mesh,&
             "ProjectedScalarViscosity")
        
        call project_field_full_Vector(velocity, continuous_velocity, positions)
        call strain_rate(continuous_velocity,positions,strain_rate_tensor)
        call tensor_second_invariant(strain_rate_tensor, strain_rate_invariant)
        
        call bound(strain_rate_invariant,0.0,huge(1.0))

        scalar_viscosity%val=&
             rheology%mu_inf+(rheology%mu_0-rheology%mu_inf)*&
             (1.0d0+rheology%lambda**2*4.0*strain_rate_invariant%val)&
             **((rheology%exponent-1.0d0)/2.0)

        call project_field_full_scalar(scalar_viscosity,&
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


      end subroutine calculate_carreau_rheology

      subroutine calculate_viscoelastic_local_rheology(positions,state,viscosity,stress,rheology,&
           continuous_velocity_mesh)
        type(tensor_field)  :: viscosity, stress
        type(viscoelastic_rheology_type) :: rheology
        type(state_type) :: state
        type(vector_field) :: positions
        type(mesh_type) :: continuous_velocity_mesh

        type(vector_field), pointer :: velocity
        type(tensor_field), pointer :: old_stress
        type(tensor_field) :: strain_rate_tensor, projected_strain_rate
        
        type(vector_field) :: continuous_velocity

        integer :: i, j, dim

        real :: dt


        velocity=>extract_vector_field(state,"Velocity")
        old_stress=> extract_tensor_field(state,"OldStress")
        dim=mesh_dim(velocity)
        call allocate(continuous_velocity,mesh=continuous_velocity_mesh,&
             name="ContinuousVelocity",dim=dim)
        call allocate(strain_rate_tensor,velocity%mesh,name="StrainRateTensor",&
             dim=[dim,dim])
        call allocate(projected_strain_rate,stress%mesh,name="ProjectedStrainRateTensor",&
             dim=[dim,dim])
        
        call project_field(velocity, continuous_velocity, positions)
        call strain_rate(continuous_velocity,positions,strain_rate_tensor)
        call project_field(strain_rate_tensor,projected_strain_rate,positions)
        
        stress%val=(old_stress%val-2.0*rheology%mu_p*projected_strain_rate%val)/(1.0+timestep/rheology%lambda)
        viscosity%val=rheology%mu_s+rheology%mu_p
        
        call deallocate(strain_rate_tensor)
        call deallocate(projected_strain_rate)
        call deallocate(continuous_velocity)

        


      end subroutine calculate_viscoelastic_local_rheology


      subroutine calculate_rheologies(state,rheology)
        type(state_type), dimension(:) :: state
        type(rheology_type), dimension(:) :: rheology

        integer :: i, j

        type(scalar_field), pointer :: cvf, pressure
        type(vector_field), pointer :: positions
        type(tensor_field), pointer :: viscosity, viscosity2, stress
        
        type(mesh_type), pointer :: continuous_velocity_mesh, material_mesh, ovmesh
        integer :: stat

        positions=>extract_vector_field(state(1),"Coordinate")
        continuous_velocity_mesh=>extract_mesh(state(1),"VelocityMesh_Continuous")
        
        if ( .not. has_mesh(state(1),"PressureMesh_Discontinuous") ) then
           pressure=> extract_scalar_field(state(1),"Pressure")
           allocate(ovmesh)
           ovmesh=make_mesh(positions%mesh,&
                shape=pressure%mesh%shape,&
              continuity=-1,name="PressureMesh_Discontinuous")
           do i=1,size(state)
              call insert(state(i),ovmesh,"PressureMesh_Discontinuous")
           end do
           call deallocate(ovmesh)
           deallocate(ovmesh)
           material_mesh=>extract_mesh(state(1),"PressureMesh_Discontinuous")
           material_mesh%option_path='/geometry/mesh::PressureMesh_Discontinuous'          
           call add_option(trim(material_mesh%option_path)// "/from_mesh",stat)
           call add_option(trim(material_mesh%option_path)// "/stat/exclude_from_stat",stat)
        else
           material_mesh=>extract_mesh(state(1),"PressureMesh_Discontinuous")
        end if
        do i=1,size(state)
           if (rheology(i)%rheology_model==2) then
              viscosity=>extract_tensor_field(state(i),"Viscosity")
              call calculate_power_law_rheology(positions,state(i),state(1),&
                   viscosity,rheology(i)%power_law_ptr,&
                   continuous_velocity_mesh,&
                   material_mesh)
           else if (rheology(i)%rheology_model==3) then
              viscosity=>extract_tensor_field(state(i),"Viscosity")
              call calculate_carreau_rheology(positions,state(i),&
                   viscosity,rheology(i)%carreau_ptr,&
                   continuous_velocity_mesh,&
                   material_mesh)
           else if (rheology(i)%rheology_model==4) then
              viscosity=>extract_tensor_field(state(i),"Viscosity")
              viscosity=>extract_tensor_field(state(i),"Stress")
              call calculate_viscoelastic_local_rheology(positions,state(i),&
                   viscosity,stress,rheology(i)%viscoelastic_ptr,continuous_velocity_mesh)
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

      subroutine update_old_full_stress(state,continuous_velocity_mesh,rheology,positions)
        type (state_type), intent(inout) :: state
        type(mesh_type), intent(in)      :: continuous_velocity_mesh
        type(viscoelastic_rheology_type)  :: rheology
        type(vector_field) :: positions

        type(tensor_field), pointer :: new_stress, old_stress, viscosity
        type(tensor_field) :: strain_rate_tensor, projected_strain_rate
        type(vector_field) :: continuous_velocity
        type(vector_field), pointer :: velocity

        integer :: dim
        
        new_stress => extract_tensor_Field(state,"Stress") 
        old_stress => extract_tensor_Field(state,"OldStress") 

        
        velocity=>extract_vector_field(state,"Velocity") 
        dim=mesh_dim(velocity)
        call allocate(continuous_velocity,mesh=continuous_velocity_mesh,&
             name="ContinuousVelocity",dim=dim)
        call allocate(strain_rate_tensor,mesh=velocity%mesh,name="StrainRateTensor",&
             dim=[dim,dim])
        call allocate(projected_strain_rate,mesh=new_stress%mesh,name="ProjectedStrainRateTensor",&
             dim=[dim,dim])
        
        call project_field(velocity, continuous_velocity, positions)
        call strain_rate(continuous_velocity,positions,strain_rate_tensor)
        call project_field(strain_rate_tensor,projected_strain_rate,positions)
        
        old_stress%val=new_stress%val+2.0*rheology%mu_p*projected_strain_rate%val
        
        call deallocate(strain_rate_tensor)
        call deallocate(projected_strain_rate)
        call deallocate(continuous_velocity)

      end subroutine update_old_full_stress

    end module multiphase_rheology
