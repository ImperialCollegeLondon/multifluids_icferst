#include "confdefs.h"

subroutine test_shear_dependent_rheology

  use fields
  use field_derivatives
  use state_module
  use spud
  use populate_state_module
  use vtk_interfaces
  use unittest_tools

  use multiphase_rheology

  implicit none

  type(state_type), pointer, dimension(:) :: state
  type(rheology_type), dimension(1) :: rheology

  type(tensor_field), pointer :: viscosity
  type(tensor_field) :: solution, error, strain_rate_tensor
  type(vector_field), pointer :: positions, velocity
  logical :: fail

  interface
     function soln(pos)
       real, dimension(:), intent(in) :: pos
       real soln(size(pos),size(pos))
     end function soln
  end interface

  call load_options("data/shear_dependent_rheology.mpml")
  call populate_state(state)
  call initialize_rheologies(state,rheology)
  call calculate_rheologies(state,rheology)

  viscosity => extract_tensor_field(state(1),"Viscosity")
  velocity => extract_vector_field(state(1),"Velocity")
  positions => extract_vector_field(state(1),"Coordinate")
  call allocate(solution,viscosity%mesh,"AnalyticSolution",dim=[3,3])
  call allocate(error,viscosity%mesh,"Error",dim=[3,3])
  call allocate(strain_rate_tensor ,viscosity%mesh,"RateOfStrain",dim=[3,3])
  

  call set_from_function(solution,soln,positions)
  
  call set(error,viscosity)
  call addto(error,solution,scale=-1.0)
  
  call vtk_write_fields("data/shear_dependent_rheology",&
       0, positions, viscosity%mesh, &
     vfields=(/ velocity /), &
     tfields=(/ viscosity,solution,error,strain_rate_tensor /))

fail = maxval( abs( error%val ))> 1e-5
  call report_test("[shear_dependent_rheology]", fail, .false., "shear dependent rheology different from analytic solution expected")

  call deallocate(solution)
  call deallocate(error)
  call deallocate(strain_rate_tensor)
  call finalize_rheologies([rheology])
  call deallocate(state)

  call report_test_no_references()

end subroutine test_shear_dependent_rheology


function soln(pos) result(mu)

  implicit none

  real, dimension(:) :: pos
  real :: mu(3,3)
  real :: x,y,z

  real :: tau_0, kappa, gamma , n

  tau_0=1.0e-3
  kappa=1.0e-2
  n=0.9
 
  x=pos(1)
  y=pos(2)
  z=pos(3)

  ! u = [y,z**2,x+y]
  
  ! D = [ [0,1,1] , [ 1 , 0, 2*z+1] , [1 ,1+2*z, 0] ]

  gamma = max(1.0e-8,sqrt((2.0+(1.0+2.0*z)**2))) 
  


  mu=tau_0/gamma+kappa*gamma**(n-1.0)

end function soln
