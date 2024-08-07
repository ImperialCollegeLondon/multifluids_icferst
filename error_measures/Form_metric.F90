#include "fdebug.h"


module form_metric_field

  use fldebug
  use vector_tools
  use global_parameters, only : OPTION_PATH_LEN
  use futils, only: int2str
  use spud
  use unittest_tools
  use metric_tools
  use fields
  use merge_tensors
  use state_module
  use vtk_interfaces, only: vtk_write_fields
  use recovery_estimator
  use field_options
  
  implicit none

  interface bound_metric
    module procedure bound_metric_anisotropic, bound_metric_minmax
  end interface

  interface form_metric
    module procedure form_metric_state, form_metric_weight
  end interface

  private
  public :: form_metric, bound_metric, p_norm_scale_metric

  contains

  subroutine form_metric_state(state, hessian, field)
    !!< Form the metric for the field field, storing it in hessian.
    type(state_type), intent(in) :: state
    type(tensor_field), intent(inout) :: hessian
    type(scalar_field), intent(inout) :: field

    type(scalar_field), pointer :: adweit

    integer :: idx, dim
    logical :: allocated

    idx = index(trim(field%name), "%")
    if (idx == 0) then
      adweit => extract_scalar_field(state, trim(field%name) // "InterpolationErrorBound")
      allocated = .false.
    else
      read(field%name(idx+1:len_trim(field%name)), *) dim
      adweit => extract_scalar_field(state, field%name(1:idx-1) // "InterpolationErrorBound%" // int2str(dim), allocated=allocated)
    end if
    
    call form_metric(hessian, field, adweit, state)

    if (allocated) then
      deallocate(adweit)
    end if

  end subroutine form_metric_state

  subroutine form_metric_weight(hessian, field, adweit, state)
    !!< Form the metric for the field field, storing it in hessian.
    type(tensor_field), intent(inout) :: hessian
    type(scalar_field), intent(inout) :: field
    type(scalar_field), intent(in) :: adweit
    type(state_type), intent(in) :: state

    integer :: p, stat

    if (have_adapt_opt(trim(field%option_path), "/adaptivity_options/relative_measure")) then
      ewrite(2,*) "Forming relative metric"
      call relative_metric(hessian, field, adweit)
    else
      ewrite(2,*) "Forming absolute metric"
            
      call get_p_norm(field, p, stat)
      if(stat == 0) then
        ewrite(2, *) "Norm degree: ", p
        call absolute_metric(hessian, adweit, p)
      else
        ewrite(2, *) "Norm degree: inf"
        call absolute_metric(hessian, adweit)
      end if
    end if

    ewrite(2,*) "Bounding metric"
    call bound_metric(hessian, state)
  
  contains
  
    subroutine get_p_norm(field, p, stat)
      type(scalar_field), intent(in) :: field
      integer, intent(out) :: p
      integer, intent(out) :: stat
      
      character(len = OPTION_PATH_LEN) :: option_path
      
      stat = 0
      
      option_path = complete_field_path(field%option_path, stat = stat)
      if(stat /= 0) return
      
      call get_option(trim(option_path) // "/adaptivity_options/absolute_measure/p_norm", p, stat = stat)
      
    end subroutine get_p_norm
  
  end subroutine form_metric_weight

  subroutine absolute_metric(metric, adweit, p)
    !!< Construct the metric using the absolute error formulation.
    type(tensor_field), intent(inout) :: metric
    type(scalar_field), intent(in) :: adweit
    !! Norm degree. Defaults to inf-norm.
    integer, optional, intent(in) :: p
    
    integer :: i
    
    if(present(p)) call p_norm_scale_metric(metric, p)
    
    do i = 1, node_count(metric)
      call set(metric, i, node_val(metric, i) / node_val(adweit, i))
    end do
      
  end subroutine absolute_metric
  
  subroutine p_norm_scale_metric(metric, p)
    !!< Apply the p-norm scaling to the metric, as in Chen Sun and Zu,
    !!< Mathematics of Computation, Volume 76, Number 257, January 2007,
    !!< pp. 179-204
    
    type(tensor_field), intent(inout) :: metric
    integer, intent(in) :: p
    
    real :: m_det
    real, dimension(mesh_dim(metric), mesh_dim(metric)) :: evecs
    real, dimension(mesh_dim(metric)) :: evals
    integer :: i, j, n
    
    assert(p > 0)
    
    n = mesh_dim(metric)
    
    do i = 1, node_count(metric)
      ! We really need det(metric) to be positive here
      ! This happens again in bound_metric, but it doesn't hurt to do it twice
      ! (Patrick assures me eigenrecompositions are cheap)
      call eigendecomposition_symmetric(node_val(metric, i), evecs, evals)
      do j = 1, n
        evals(j) = abs(evals(j))
      end do
      call eigenrecomposition(metric%val(:, :, i), evecs, evals)
    
      m_det = 1.0
      do j = 1, size(evals)
        m_det = m_det * evals(j)
      end do
      m_det = max(m_det, epsilon(0.0))
      
      call set(metric, i, node_val(metric, i) * (m_det ** (-1.0 / (2.0 * p + n))))
    end do
    
  end subroutine p_norm_scale_metric

  subroutine relative_metric(metric, field, adweit, p)
    !!< Construct the metric using the relative error formulation.
    type(tensor_field), intent(inout) :: metric
    type(scalar_field), intent(in) :: field, adweit
    real :: minfield, maxfield
    !! Norm degree. Defaults to inf-norm.
    integer, optional, intent(in) :: p

    integer :: i
    real, parameter :: tolerance = 1e-8

    if(present(p)) call p_norm_scale_metric(metric, p)

    call field_stats(field, min=minfield)
    call field_stats(field, max=maxfield)

    do i = 1, node_count(metric)
      call set(metric, i, node_val(metric, i) / ( node_val(adweit, i)*(maxfield-minfield) + tolerance ))
    end do

  end subroutine relative_metric

  subroutine bound_metric_anisotropic(hessian, state, stat)
    !!< Implement the anisotropic edge length bounds.
    type(tensor_field), intent(inout) :: hessian
    type(state_type), intent(in) :: state
    integer, optional :: stat

    type(tensor_field), pointer :: min_bound, max_bound

    ! min and max refer to edge lengths, not eigenvalues

    if(present(stat)) stat = 0

    min_bound => extract_tensor_field(state, "MaxMetricEigenbound", stat=stat)
    max_bound => extract_tensor_field(state, "MinMetricEigenbound", stat=stat)
    if(present(stat)) then
      if(stat/=0) return
    end if
    if(any(min_bound%dim /= hessian%dim).or.any(max_bound%dim /= hessian%dim)) then
      if(present(stat)) then
        stat = 1
        return
      else
        FLAbort("Incompatible tensor dimensions")
      end if
    end if

    call bound_metric(hessian, min_bound, max_bound)

  end subroutine bound_metric_anisotropic
  
  subroutine bound_metric_minmax(hessian, min_bound, max_bound)
    !!< Implement the anisotropic edge length bounds.
    type(tensor_field), intent(inout) :: hessian
    type(tensor_field), intent(in) :: min_bound, max_bound

    real, dimension(hessian%dim(1), hessian%dim(2)) :: max_tensor, min_tensor, evecs
    real, dimension(hessian%dim(1)) :: evals
    integer :: i, j

    ! min and max refer to edge lengths, not eigenvalues
    do i=1,node_count(hessian)
      call eigendecomposition_symmetric(hessian%val(:, :, i), evecs, evals)
      do j=1,hessian%dim(1)
        evals(j) = abs(evals(j))
      end do
      call eigenrecomposition(hessian%val(:, :, i), evecs, evals)
      assert(all(evals >= 0.0))

      max_tensor = node_val(max_bound, i)
      min_tensor = node_val(min_bound, i)
      call merge_tensor(hessian%val(:, :, i), max_tensor)
      call merge_tensor(hessian%val(:, :, i), min_tensor, aniso_min=.true.)
    end do

  end subroutine bound_metric_minmax

end module form_metric_field
