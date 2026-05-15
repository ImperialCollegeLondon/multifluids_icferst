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

module python_diagnostics

  use fldebug
  use global_parameters, only : PYTHON_FUNC_LEN, OPTION_PATH_LEN
  use spud
  use fields
  use python_state
  use state_module
  use parallel_tools, only: getprocno

  implicit none
  
  private
  
  public :: calculate_scalar_python_diagnostic, &
    & calculate_scalar_python_diagnostic_refreshed, &
    & calculate_vector_python_diagnostic, calculate_tensor_python_diagnostic
    
contains

  subroutine calculate_scalar_python_diagnostic(states, state_index, s_field, current_time, dt)
    !!< Set a field from Python
    !!< So add the whole state and make a variable with the diagnostic
    !!< field available to the interpreter
  
    type(state_type), dimension(:), target, intent(inout) :: states
    integer, intent(in) :: state_index
    type(scalar_field), intent(inout) :: s_field
    real, intent(in) :: current_time
    real, intent(in) :: dt

#ifdef HAVE_NUMPY    
    character(len = PYTHON_FUNC_LEN) :: pycode
    character(len = 30) :: buffer
    character(len = OPTION_PATH_LEN) :: material_phase_support
    type(state_type), pointer :: this_state
    
    ewrite(2,*) 'in calculate_scalar_python_diagnostic'
    ! Clean up to make sure that nothing else interferes
    call python_reset()
    
    call get_option(trim(s_field%option_path)&
         //"/diagnostic/algorithm/material_phase_support",material_phase_support)

    ewrite(2,*) 'material_phase_support: '//trim(material_phase_support)
    select case(material_phase_support)
    case("single")
       call python_add_state(states(state_index))

    case ("multiple")
       call python_add_states(states)
       this_state=>states(state_index)
       call python_run_string("state = states['"//trim(this_state%name)//"']")

    case default
       ewrite(0,*) trim(material_phase_support)&
            //" is not a valid value for material_phase_support"
       FLExit("Options file error")
    end select

    call python_run_string("field = state.scalar_fields['"//trim(s_field%name)//"']")
    write(buffer,*) current_time
    call python_run_string("time="//trim(buffer))
    write(buffer,*) dt
    call python_run_string("dt="//trim(buffer))  
      
    ! And finally run the user's codey
    call get_option(trim(s_field%option_path)//"/diagnostic/algorithm",pycode)
    call python_run_string(trim(pycode))
    
    ! Cleanup
    call python_reset()
#else
    FLAbort("Python diagnostic fields require NumPy, which cannot be located.")
#endif

    ewrite(2,*) 'leaving calculate_scalar_python_diagnostic'

  end subroutine calculate_scalar_python_diagnostic

  subroutine calculate_scalar_python_diagnostic_refreshed(states, state_index, s_field, current_time, dt)
    type(state_type), dimension(:), target, intent(inout) :: states
    integer, intent(in) :: state_index
    type(scalar_field), intent(inout) :: s_field
    real, intent(in) :: current_time
    real, intent(in) :: dt

#ifdef HAVE_NUMPY
    character(len = PYTHON_FUNC_LEN) :: pycode
    character(len = 30) :: buffer
    character(len = 10) :: nbuf, rankbuf
    type(scalar_field), pointer :: sf
    type(vector_field), pointer :: vf
    integer :: i, j, ndim, nnodes, nfields, io_unit, myrank
    character(len=256) :: tmpfile, resfile

    io_unit = 4321
    myrank = getprocno() - 1
    write(rankbuf, '(I0)') myrank

    ! Define mock classes
    call python_run_string( &
        'import numpy as np' // ACHAR(10) // &
        'class _SF:' // ACHAR(10) // &
        '  def __init__(s,d):' // ACHAR(10) // &
        '    s._d=np.array(d,dtype="float64");s.node_count=len(d)' // ACHAR(10) // &
        '  def node_val(s,i): return float(s._d[i])' // ACHAR(10) // &
        '  def set(s,i,v): s._d[i]=v' // ACHAR(10) // &
        'class _VF:' // ACHAR(10) // &
        '  def __init__(s,comps):' // ACHAR(10) // &
        '    s._c=[np.array(c,dtype="float64") for c in comps]' // ACHAR(10) // &
        '    s.node_count=len(s._c[0])' // ACHAR(10) // &
        '  def node_val(s,i): return [float(c[i]) for c in s._c]' // ACHAR(10) // &
        'class _ST:' // ACHAR(10) // &
        '  def __init__(s): s.scalar_fields={};s.vector_fields={}' // ACHAR(10) // &
        'state=_ST()')

    ! Scalar fields via binary file (rank-specific)
    nfields = scalar_field_count(states(state_index))
    do i = 1, nfields
        sf => extract_scalar_field(states(state_index), i)
        nnodes = size(sf%val)
        write(nbuf, '(I0)') nnodes
        tmpfile = '/tmp/icferst_sf_' // trim(rankbuf) // '.bin'
        open(unit=io_unit, file=trim(tmpfile), form='unformatted', &
             access='stream', status='replace', action='write')
        write(io_unit) sf%val(1:nnodes)
        close(io_unit)
        call python_run_string( &
            "state.scalar_fields['" // trim(sf%name) // &
            "']=_SF(np.fromfile('" // trim(tmpfile) // &
            "',dtype='float64',count=" // trim(nbuf) // "))")
    end do

    ! Vector fields via binary file (rank-specific)
    nfields = vector_field_count(states(state_index))
    do i = 1, nfields
        vf => extract_vector_field(states(state_index), i)
        ndim = vf%dim
        nnodes = vf%mesh%nodes
        write(nbuf, '(I0)') nnodes
        do j = 1, ndim
            tmpfile = '/tmp/icferst_vf_' // trim(rankbuf) // '.bin'
            open(unit=io_unit, file=trim(tmpfile), form='unformatted', &
                 access='stream', status='replace', action='write')
            write(io_unit) vf%val(j, 1:nnodes)
            close(io_unit)
            write(buffer, '(I0)') j
            call python_run_string( &
                "_vc" // trim(buffer) // &
                "=np.fromfile('" // trim(tmpfile) // &
                "',dtype='float64',count=" // trim(nbuf) // ")")
        end do
        if (ndim == 1) then
            call python_run_string("state.vector_fields['" // &
                trim(vf%name) // "']=_VF([_vc1])")
        else if (ndim == 2) then
            call python_run_string("state.vector_fields['" // &
                trim(vf%name) // "']=_VF([_vc1,_vc2])")
        else
            call python_run_string("state.vector_fields['" // &
                trim(vf%name) // "']=_VF([_vc1,_vc2,_vc3])")
        end if
    end do

    ! Output field
    nnodes = size(s_field%val)
    write(nbuf, '(I0)') nnodes
    tmpfile = '/tmp/icferst_sf_' // trim(rankbuf) // '.bin'
    open(unit=io_unit, file=trim(tmpfile), form='unformatted', &
         access='stream', status='replace', action='write')
    write(io_unit) s_field%val(1:nnodes)
    close(io_unit)
    call python_run_string( &
        "field=_SF(np.fromfile('" // trim(tmpfile) // &
        "',dtype='float64',count=" // trim(nbuf) // "))")

    ! Time and dt
    write(buffer, *) current_time
    call python_run_string("time=" // trim(buffer))
    write(buffer, *) dt
    call python_run_string("dt=" // trim(buffer))

    ! Step 6: Run user code
    call get_option(trim(s_field%option_path) &
         // "/diagnostic/algorithm", pycode)
    call python_run_string(trim(pycode))

    ! Write results back via binary file
    resfile = '/tmp/icferst_res_' // trim(rankbuf) // '.bin'
    call python_run_string( &
        "field._d.astype('float64').tofile('" // trim(resfile) // "')")
    open(unit=io_unit, file=trim(resfile), form='unformatted', &
         access='stream', status='old', action='read')
    read(io_unit) s_field%val(1:nnodes)
    close(io_unit)

#else
    FLAbort("Requires NumPy.")
#endif

  end subroutine calculate_scalar_python_diagnostic_refreshed

  subroutine calculate_vector_python_diagnostic(states, state_index, v_field, current_time, dt)
    !!< Set a field from Python
    !!< So add the whole state and make a variable with the diagnostic
    !!< field available to the interpreter
    
    type(state_type), dimension(:), target, intent(inout) :: states
    integer, intent(in) :: state_index
    type(vector_field), intent(inout) :: v_field
    real, intent(in) :: current_time
    real, intent(in) :: dt
    
#ifdef HAVE_NUMPY
    character(len = PYTHON_FUNC_LEN) :: pycode
    character(len = 30) :: buffer
    character(len = OPTION_PATH_LEN) :: material_phase_support
    type(state_type), pointer :: this_state
    
    ! Clean up to make sure that nothing else interferes
    call python_reset()
    
    call get_option(trim(v_field%option_path)&
         //"/diagnostic/algorithm/material_phase_support",material_phase_support)

    select case(material_phase_support)
    case("single")
       call python_add_state(states(state_index))

    case ("multiple")
       call python_add_states(states)
       this_state=>states(state_index)
       call python_run_string("state = states['"//trim(this_state%name)//"']")

    case default
       ewrite(0,*) trim(material_phase_support)&
            //" is not a valid value for material_phase_support"
       FLExit("Options file error")
    end select

    call python_run_string("field = state.vector_fields['"//trim(v_field%name)//"']")
    write(buffer,*) current_time
    call python_run_string("time="//trim(buffer))
    write(buffer,*) dt
    call python_run_string("dt="//trim(buffer))  
      
    ! And finally run the user's code
    call get_option(trim(v_field%option_path)//"/diagnostic/algorithm",pycode)
    call python_run_string(trim(pycode))
    
    ! Cleanup
    call python_reset()
#else
    FLAbort("Python diagnostic fields require NumPy, which cannot be located.")
#endif
    
  end subroutine calculate_vector_python_diagnostic
  
  subroutine calculate_tensor_python_diagnostic(states, state_index, t_field, current_time, dt)
    !!< Set a field from Python
    !!< So add the whole state and make a variable with the diagnostic
    !!< field available to the interpreter
    
    type(state_type), dimension(:), target, intent(inout) :: states
    integer, intent(in) :: state_index
    type(tensor_field), intent(inout) :: t_field
    real, intent(in) :: current_time
    real, intent(in) :: dt
    
#ifdef HAVE_NUMPY
    character(len = PYTHON_FUNC_LEN) :: pycode
    character(len = 30) :: buffer
    character(len = OPTION_PATH_LEN) :: material_phase_support
    type(state_type), pointer :: this_state
    
    ! Clean up to make sure that nothing else interferes
    call python_reset()
    
    call get_option(trim(t_field%option_path)&
         //"/diagnostic/algorithm/material_phase_support",material_phase_support)

    select case(material_phase_support)
    case("single")
       call python_add_state(states(state_index))

    case ("multiple")
       call python_add_states(states)
       this_state=>states(state_index)
       call python_run_string("state = states['"//trim(this_state%name)//"']")

    case default
       ewrite(0,*) trim(material_phase_support)&
            //" is not a valid value for material_phase_support"
       FLExit("Options file error")
    end select

    call python_run_string("field = state.tensor_fields['"//trim(t_field%name)//"']")
    write(buffer,*) current_time
    call python_run_string("time="//trim(buffer))
    write(buffer,*) dt
    call python_run_string("dt="//trim(buffer))  
      
    ! And finally run the user's code
    call get_option(trim(t_field%option_path)//"/diagnostic/algorithm",pycode)
    call python_run_string(trim(pycode))
    
    ! Cleanup
    call python_reset()
#else
    FLAbort("Python diagnostic fields require NumPy, which cannot be located.")
#endif
    
  end subroutine calculate_tensor_python_diagnostic

end module python_diagnostics
