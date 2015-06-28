#include "fdebug.h"

module limit_metric_module


  use elements
  use fields
  use data_structures
  use fldebug
  use meshdiagnostics
  use spud
  use parallel_fields
  use vector_tools, only : determinant => det

  implicit none

  logical :: limit_by_region=.false.

  private

  public :: limit_metric, limit_metric_elements, expected_elements, &
    & expected_nodes, determinant, limit_by_region, build_region_list

  interface expected_elements
    module procedure expected_elements_whole_mesh, expected_elements_by_region
  end interface expected_elements

  interface expected_nodes
    module procedure expected_nodes_expected_elements, expected_nodes_expected_elements_by_region, expected_nodes_metric
  end interface expected_nodes
  
  interface limit_metric
    module procedure limit_metric_nodes_options, limit_metric_nodes_minmax, &
      & limit_metric_nodes_minmax_by_region, limit_metric_nodes_target
  end interface limit_metric
  
  interface limit_metric_elements
    module procedure limit_metric_elements_minmax, limit_metric_elements_target
  end interface limit_metric_elements

contains

  subroutine limit_metric_nodes_options(positions, metric)
    type(tensor_field), intent(inout) :: metric
    type(vector_field), intent(in) :: positions

    character(len = *), parameter :: base_path = "/mesh_adaptivity/hr_adaptivity"
    integer :: nodes, stat
    real :: increase_tolerance

    type(integer_set), dimension(:), allocatable :: region_list
    integer, dimension(:), allocatable :: max_nodes, min_nodes
    integer :: no_of_regions

    call build_region_list(region_list, max_nodes, min_nodes)
    no_of_regions=size(region_list)
    
    call mesh_stats(positions, nodes = nodes)

    call get_option(base_path // "/minimum_number_of_nodes", min_nodes(no_of_regions+1), default = 1)
    if(have_option(base_path // "/minimum_number_of_nodes/per_process")) then
      min_nodes(no_of_regions+1) = min_nodes(no_of_regions+1) * getnprocs()
    end if
    call get_option(base_path // "/maximum_number_of_nodes", max_nodes(no_of_regions+1), default = 100000)
    if(have_option(base_path // "/maximum_number_of_nodes/per_process")) then
      max_nodes(no_of_regions+1) = max_nodes(no_of_regions+1) * getnprocs()
    end if
    call get_option(base_path // "/max_node_increase", increase_tolerance, stat = stat)
    if(stat == SPUD_NO_ERROR) then
      max_nodes(no_of_regions) = min(max_nodes(no_of_regions+1), int(nodes * increase_tolerance))
    end if

    if (no_of_regions>0) then
       max_nodes(no_of_regions+1)=max_nodes(no_of_regions+1)-sum(max_nodes(1:no_of_regions))
       min_nodes(no_of_regions+1)=min_nodes(no_of_regions+1)-sum(min_nodes(1:no_of_regions))
       if (min_nodes(no_of_regions+1)<1) min_nodes(no_of_regions+1)=1
       limit_by_region=.true.
    else
       limit_by_region=.false.
    end if


    call limit_metric(positions, metric, region_list,min_nodes, max_nodes)

    call deallocate(region_list)

  end subroutine limit_metric_nodes_options

  subroutine limit_metric_nodes_minmax(positions, metric, min_nodes, max_nodes)
    type(vector_field), intent(in) :: positions
    type(tensor_field), intent(inout) :: metric
    integer, intent(in):: min_nodes
    integer, intent(in) :: max_nodes


    type(integer_set), dimension(0) :: region_list


    call limit_metric_nodes_minmax_by_region(positions, metric, region_list, [min_nodes], [max_nodes])

  end subroutine limit_metric_nodes_minmax
  
  subroutine limit_metric_nodes_minmax_by_region(positions, metric, region_list, min_nodes, max_nodes)
    type(vector_field), intent(in) :: positions
    type(tensor_field), intent(inout) :: metric
    type(integer_set), dimension(:) :: region_list
    integer, intent(in), dimension(:) :: min_nodes
    integer, intent(in), dimension(:) :: max_nodes

    integer, dimension(size(min_nodes)) :: elements, max_eles, min_eles, nodes
    ! The ratio of elements to nodes
    real, dimension(size(min_nodes)) :: eles_per_node


    integer, dimension(:,:), allocatable :: work
    integer :: ele, region, i
    logical :: found_region
    integer,  pointer, dimension(:) :: nod
    
    
    assert(all(min_nodes > 0))
    assert(sum(max_nodes) >= sum(min_nodes))
    
    if (size(region_list)==0) then
       call mesh_stats(positions, nodes = nodes(1), elements = elements(1))
    else
       allocate(work(node_count(positions),size(nodes)))
       work=0
       elements=0
       do ele=1,ele_count(positions)
          found_region=.false.
          nod=>ele_nodes(positions,ele)
          do region=1,size(region_list)
             if (has_value(region_list(region),&
                  ele_region_id(positions,ele))) then
                do i=1,size(nod)
                   if (node_owned(positions,nod(i))) &
                        work(nod(i),region)=1
                end do
                if (element_owned(positions,ele)) &
                     elements(region)=elements(region)+1
                found_region=.true.
                exit
             end if
          end do
          if (.not. found_region) then
             region=size(elements)
             do i=1,size(nod)
                if (node_owned(positions,nod(i))) &
                     work(nod(i),region)=1
             end do
             if (element_owned(positions,ele)) &
                  elements(region)=elements(region)+1
          end if
       end do
       nodes=sum(work,dim=1)
       call allsum(nodes)
       call allsum(elements)
    end if
    
    ! FIXME: maybe a better way to do this?
    eles_per_node = float(elements) / float(nodes)
    
    min_eles = int(eles_per_node * min_nodes)
    max_eles = int(eles_per_node * max_nodes)

    if(any(min_eles < 0)) then
      ewrite(-1, *) "Minimum elements: ", min_eles
      FLAbort("Invalid minimum number of elements")
    end if
    if(any(max_eles < 0)) then
      ewrite(-1, *) "Maximum elements: ", max_eles
      FLAbort("Invalid maximum number of elements")
    end if

    call limit_metric_elements(positions, metric, region_list,min_eles, max_eles)
    
  end subroutine limit_metric_nodes_minmax_by_region
  
  subroutine limit_metric_elements_minmax(positions, metric,region_list, min_eles, max_eles)
    type(vector_field), intent(in) :: positions
    type(tensor_field), intent(inout) :: metric
    type(integer_set), dimension(:) :: region_list
    integer, intent(in), dimension(:) :: min_eles
    integer, intent(in), dimension(:) :: max_eles
    
    integer, dimension(size(min_eles)) :: expected_eles
    ! the scaling factor to divide the metric by
    real, dimension(size(min_eles)) :: beta
    type(scalar_field) :: beta_field
    integer :: i,j, ele, region
    real :: val
    integer, pointer, dimension(:) :: nod
    logical :: found_region 
    
    expected_eles = expected_elements(positions, metric, region_list, global = .true.)
    
    do j=1,size(beta)
       if (expected_eles(j) > max_eles(j)) then
          beta = ((1.0 / expected_eles(j)) * max_eles(j)) ** (2.0 / positions%dim)
          ewrite(2,*) "Scaling factor to conform to maximum node limit: ", beta(j)
       else if(expected_eles(j) < min_eles(j)) then
          beta = ((1.0 / expected_eles(j)) * min_eles(j)) ** (2.0 / positions%dim)
          ewrite(2,*) "Scaling factor to conform to minimum node limit: ", beta(j)
       else
          beta(j)=1.0
       end if
    end do

    if (size(beta)==1 .and. beta(1) .ne. 1.0) then
       call scale(metric, beta(1))
    else if (size(beta)>1) then
       call allocate(beta_field,metric%mesh,"ScaleFactor")
       call zero(beta_field)
       
       do ele=1, ele_count(metric)

          found_region=.false.
          do region=1,size(region_list)
             if (has_value(region_list(region),ele_region_id(metric,ele))) then
                nod=>ele_nodes(metric,ele)
                do i=1, size(nod)
                   val=node_val(beta_field,nod(i))
                   if (val<beta(region)) &
                        call set(beta_field,nod(i),beta(region))
                end do
                found_region=.true.
             end if
          end do
          if (.not. found_region) then
             region=size(region_list)+1
             nod=>ele_nodes(metric,ele)
             do i=1, size(nod)
                val=node_val(beta_field,nod(i))
                if (val<beta(region)) &
                     call set(beta_field,nod(i),beta(region))
             end do
          end if
       end do            

       call scale(metric,beta_field)
       call deallocate(beta_field)

    end if
    
  end subroutine limit_metric_elements_minmax
  
  subroutine limit_metric_nodes_target(positions, metric, region_list, target_nodes)
    type(vector_field), intent(in) :: positions
    type(tensor_field), intent(inout) :: metric
    type(integer_set), dimension(:) :: region_list
    integer, intent(in), dimension(:) :: target_nodes
    
    call limit_metric(positions, metric, region_list,min_nodes = target_nodes, max_nodes = target_nodes)
    
  end subroutine limit_metric_nodes_target
  
  subroutine limit_metric_elements_target(positions, metric, region_list,target_eles)
    type(vector_field), intent(in) :: positions
    type(tensor_field), intent(inout) :: metric
    type(integer_set), dimension(:) :: region_list
    integer, intent(in), dimension(:) :: target_eles
    
    call limit_metric_elements(positions, metric, region_list, min_eles = target_eles, max_eles = target_eles)
    
  end subroutine limit_metric_elements_target

  subroutine build_region_list(region_list,max_nodes, min_nodes) 
  
    type(integer_set), dimension(:), allocatable, intent(inout) :: region_list
    integer, dimension(:), allocatable, intent(inout) :: max_nodes, min_nodes

    character(len = *), parameter :: base_path = "/mesh_adaptivity/hr_adaptivity"
    integer :: no_of_regions, region, id, id_shape(2)
    integer, dimension(:), allocatable :: region_ids

    no_of_regions=option_count(base_path // "/number_of_nodes_by_mesh_region" )

    allocate(region_list(no_of_regions),max_nodes(no_of_regions+1),min_nodes(no_of_regions+1))
    call allocate(region_list)
    do region=1, no_of_regions
       id_shape=option_shape(base_path // "/number_of_nodes_by_mesh_region/region_ids")
       allocate(region_ids(id_shape(1)))
       call get_option(base_path // "/number_of_nodes_by_mesh_region["//int2str(region-1)//"]/region_ids",&
            region_ids)
       do id=1,size(region_ids)
          call insert(region_list(id),region_ids(id))
       end do
       call get_option(base_path // "/number_of_nodes_by_mesh_region["//int2str(region-1)//"]/maximum_number_of_nodes",max_nodes(region))
       call get_option(base_path // "/number_of_nodes_by_mesh_region["//int2str(region-1)//"]/minimum_number_of_nodes",min_nodes(region),default=1)
       deallocate(region_ids)
       if (min_nodes(region)>max_nodes(region)) then
          FLAbort("Minimum number of nodes exceeds maximum number for at least one mesh region. Please fix your adaptivity options")
       end if

    end do

  end subroutine build_region_list

  function expected_elements_whole_mesh(old_positions, metric,  global) result(xpct)
    type(vector_field), intent(in) :: old_positions
    type(tensor_field), intent(in) :: metric
    !! If present and .true., calculate the global number of expected elements
    logical, optional, intent(in) :: global

    integer xpct, xpct_vec(1)

    type(integer_set) :: region_list(0)

    xpct_vec=expected_elements_by_region(old_positions, metric, region_list, global) 
    xpct=xpct_vec(1)

  end function expected_elements_whole_mesh

  function expected_elements_by_region(old_positions, metric, region_list, global) result(xpct)
    type(vector_field), intent(in) :: old_positions
    type(tensor_field), intent(in) :: metric
    !! If present and .true., calculate the global number of expected elements
    logical, optional, intent(in) :: global
    type(integer_set), dimension(:) :: region_list
    integer :: ele
    real, dimension(mesh_dim(metric), mesh_dim(metric)) :: avg_metric
    real :: sumvol(size(region_list)+1), det
    integer :: xpct(size(region_list)+1), region
    real :: gamma
    
    logical :: lglobal, found_region
    integer :: no_of_regions

    no_of_regions=size(region_list)

    lglobal = present_and_true(global)

    sumvol = 0.0

    ! Gamma is the volume of an optimal element
    ! (in metric space)
    select case(mesh_dim(metric))
      case(3)
        gamma = 1.0 / sqrt(72.0)
      case(2)
        gamma = sqrt(3.0) / 4.0
      case(1)
        gamma = 1.0
      case default
        FLAbort("Invalid dimension")
    end select

    do ele=1,ele_count(old_positions)
      if(lglobal) then
        if(.not. element_owned(old_positions, ele)) cycle
      end if
      avg_metric = sum(ele_val(metric, ele), 3) / ele_loc(metric, ele)
      det = determinant(avg_metric)
      found_region=.false.
      do region=1,no_of_regions
         if (has_value(region_list(region),ele_region_id(old_positions,ele))) then
            sumvol(region) = sumvol(region) &
                 + abs(sqrt(det) * simplex_volume(old_positions, ele))
            found_region=.true.
            exit
         end if
      end do
      if (.not. found_region) &
           sumvol(size(sumvol)) = sumvol(size(sumvol)) &
                 + abs(sqrt(det) * simplex_volume(old_positions, ele))
    end do
    if(lglobal) call allsum(sumvol)

    if(any(sumvol < 0.0)) then
      ewrite(-1, *) "Total volume in metric spaces: ", sumvol
      FLAbort("Negative volume")
    end if
    if(any((sumvol/gamma)>=huge(xpct(1)))) then
      ewrite(-1, *) "ERROR: The error metric &
        & indicates that number of elements required is ", sumvol/gamma
      ewrite(-1, *) "If this is what you want then, great, congratulations, &
        &this is a record. Please ask the developers to get rid of this &
        &32bit integer that's causing trouble. Otherwise, please review your &
        &error targets"
        FLExit("integer overflow")
    end if
    xpct = int(sumvol / gamma)
    if (sum(xpct) == 0) xpct(size(xpct)) = 1

    if(any(xpct < 0)) then
      ewrite(-1, *) "Expected elements: ", xpct
      FLAbort("Invalid number of expected elements")
    end if
    if (lglobal) then
      ewrite(2, *) "Expected global n/o elements: ", xpct
    else
      ewrite(2, *) "Expected n/o elements: ", xpct
    end if

  end function expected_elements_by_region

  function expected_nodes_expected_elements(old_positions, expected_eles, global) result(expected_nods)
  !!< Return the expected number of nodes based upon the supplied expected
    !!< number of elements

    type(vector_field), intent(in) :: old_positions
    integer, intent(in)  :: expected_eles
    !! If present and .true., calculate the global number of expected elemen
    logical, optional, intent(in) :: global


    type(integer_set), dimension(0) :: region_list

    integer :: expected_nods, expected_nodes_vec(1)

    expected_nodes_vec=expected_nodes_expected_elements_by_region(old_positions, [expected_eles],region_list, global)
    expected_nods=expected_nodes_vec(1)

  end function expected_nodes_expected_elements

  function expected_nodes_expected_elements_by_region(old_positions, expected_eles,region_list, global) result(expected_nods)
    !!< Return the expected number of nodes based upon the supplied expected
    !!< number of elements

    type(vector_field), intent(in) :: old_positions
    integer, intent(in), dimension(:) :: expected_eles
    !! If present and .true., calculate the global number of expected elements
    type(integer_set), dimension(:) :: region_list
    logical, optional, intent(in) :: global

    integer, dimension(size(expected_eles)) :: expected_nods, elements, nodes
    
    integer, dimension(:,:), allocatable :: work
    integer :: ele, region, i
    logical :: found_region
    integer,  pointer, dimension(:) :: nod

    if (size(region_list)==0) then
       if(present_and_true(global)) then
          call mesh_stats(old_positions, nodes = nodes(1), elements = elements(1))
       else
          nodes(1) = node_count(old_positions)
          elements(1) = ele_count(old_positions)
       end if
    else
       allocate(work(node_count(old_positions),size(nodes)))
       work=0
       elements=0
       do ele=1,ele_count(old_positions)
          found_region=.false.
          nod=>ele_nodes(old_positions,ele)
          do region=1,size(region_list)
             if (has_value(region_list(region),&
                  ele_region_id(old_positions,ele))) then
                do i=1,size(nod)
                   if (node_owned(old_positions,nod(i))) &
                        work(nod(i),region)=1
                end do
                if (element_owned(old_positions,ele)) &
                     elements(region)=elements(region)+1
                found_region=.true.
                exit
             end if
          end do
          if (.not. found_region) then
             region=size(elements)
             do i=1,size(nod)
                if (node_owned(old_positions,nod(i))) &
                     work(nod(i),region)=1
             end do
             if (element_owned(old_positions,ele)) &
                  elements(region)=elements(region)+1
          end if
       end do
       nodes=sum(work,dim=1)
       call allsum(nodes)
       call allsum(elements)
    end if
      
    ! FIXME: maybe a better way to do this?
    expected_nods = int((float(nodes) / float(elements)) * expected_eles)

  end function expected_nodes_expected_elements_by_region

  function expected_nodes_metric(old_positions, metric, region_list, global) result(expected_nods)
    !!< Return the expected number of nodes based upon the supplied metric

    type(vector_field), intent(in) :: old_positions
    type(tensor_field), intent(in) :: metric
    type(integer_set), dimension(:) :: region_list
    !! If present and .true., calculate the global number of expected elements
    logical, optional, intent(in) :: global

    integer :: expected_nods(size(region_list)+1)

    integer :: expected_eles(size(region_list)+1)

    expected_eles = expected_elements(old_positions, metric, region_list, global)
    expected_nods = expected_nodes(old_positions, expected_eles, region_list, global)

  end function expected_nodes_metric
  
  subroutine limit_metric_module_check_options()
    !!< Check metric limiting specific options
    
    character(len = *), parameter :: base_path = "/mesh_adaptivity/hr_adaptivity"
    integer :: max_nodes, min_nodes, stat
   
    if(.not. have_option(base_path)) then
      ! Nothing to check
      return
    end if
   
    call get_option(base_path // "/minimum_number_of_nodes", min_nodes, default = 0)
    call get_option(base_path // "/maximum_number_of_nodes", max_nodes, stat = stat)
    if(stat /= SPUD_NO_ERROR) then
      FLExit("Maximum number of nodes must be specified when using hr mesh adaptivity")
    else if(min_nodes > max_nodes) then
      FLExit("The minimum number of nodes cannot be greater than the maximum number of nodes")
    end if
    
  end subroutine limit_metric_module_check_options

end module limit_metric_module
