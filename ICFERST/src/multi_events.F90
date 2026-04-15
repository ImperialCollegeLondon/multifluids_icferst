!    Copyright (C) 2020 Imperial College London and others.
!
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Affero General Public License
!    as published by the Free Software Foundation,
!    version 3.0 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA
#include "fdebug.h"
!----------------------------------------------------------------------------------------
!> @brief Simulation event scheduling module for IC-FERST.
!>
!> Reads a CSV file of event times and ensures the timestep lands exactly on each event.
!>
!> Activated by adding a <simulation_events> block under <timestepping> Diamond.
!> Example Diamond:
!>
!>   <simulation_events>
!>     <from_file>
!>       <string_value>schedule.csv</string_value>
!>     </from_file>
!>     <period>
!>       <real_value rank="0">31536000</real_value>
!>     </period>
!>     <post_event_timestep>
!>       <real_value rank="0">86400</real_value>
!>     </post_event_timestep>
!>   </simulation_events>
!>
!> The general maximum timestep should be set via the adaptive timestepping
!> options in Diamond (e.g. adaptive_timestep/maximum_timestep or
!> adaptive_timestep_nonlinear/max_timestep). This module only reduces dt
!> further for event snapping and the post-event cap.
!>
!> The CSV file must have a single column with a header row, e.g.:
!>
!>   time_seconds
!>   86400
!>   172800
!>   604800
!>
!> @author Meissam Bahlali
!----------------------------------------------------------------------------------------
module multi_events

    use spud
    use fldebug
    use parallel_tools

    implicit none
    private

    public :: load_simulation_events
    public :: apply_event_timestep

    
    real, allocatable, dimension(:), save :: event_times_s ! event times in seconds (absolute or within one period if periodic)
    integer, save :: n_events = 0 ! number of events loaded from csv
    logical, save :: events_active = .false. ! whether the events module is active (csv successfully loaded)
    real, save :: event_period_s = -1.0 ! period of event repetition in seconds. Negative means non-periodic
    real, parameter :: EVENT_TOL = 1.0e-6 ! tolerance for "are we exactly on an event?"
    real, save    :: post_event_dt = -1.0 ! Timestep cap to apply immediately after landing on an event. If negative = disabled.
    logical, save :: just_landed_on_event = .false. ! Flag set to true when we just snapped to an event so next call can cap dt

contains

    !> @brief Load simulation events from the CSV file specified in the Diamond input.
    !>
    !> Must be called once before the time loop. If the <simulation_events> block is
    !> absent in Diamond, this subroutine returns immediately without doing anything.
    subroutine load_simulation_events()
        implicit none

        character(len=256) :: csv_file
        integer            :: io_unit, io_stat
        integer            :: n_lines, i
        real               :: t_val
        character(len=512) :: header_line

        ! Nothing to do if the user has not asked for event scheduling
        if (.not. have_option('/timestepping/simulation_events')) return

        call get_option('/timestepping/simulation_events/from_file', csv_file)
        call get_option('/timestepping/simulation_events/period', event_period_s, default = -1.0)
        call get_option('/timestepping/simulation_events/post_event_timestep', &
            post_event_dt, default=-1.0)

        io_unit = 99
        open(unit=io_unit, file=trim(adjustl(csv_file)), status='old', &
             action='read', iostat=io_stat)
        if (io_stat /= 0) then
            FLExit("simulation_events: Cannot open schedule file: " // trim(csv_file))
        end if

        read(io_unit, '(A)', iostat=io_stat) header_line  ! skip header
        if (io_stat /= 0) then
            close(io_unit)
            FLExit("simulation_events: Schedule file is empty (no header): " // trim(csv_file))
        end if

        n_lines = 0
        do
            read(io_unit, *, iostat=io_stat) t_val
            if (io_stat /= 0) exit
            n_lines = n_lines + 1
        end do
        close(io_unit)

        if (n_lines == 0) then
            if (getprocno() == 1) then
                ewrite(0,*) "WARNING simulation_events: Schedule file has no event rows: ", trim(csv_file)
            end if
            return
        end if

        ! Allocate and fill
        allocate(event_times_s(n_lines))
        n_events = n_lines

        open(unit=io_unit, file=trim(adjustl(csv_file)), status='old', &
             action='read', iostat=io_stat)
        read(io_unit, '(A)') header_line  ! skip header
        do i = 1, n_events
            read(io_unit, *) event_times_s(i)
        end do
        close(io_unit)

        events_active = .true.

        ! Report to  user
        if (getprocno() == 1) then
            ewrite(0,*)
            ewrite(0,*) "=================================================="
            ewrite(0,*) " simulation_events: Event scheduling activated"
            ewrite(0,*) "   File    : ", trim(csv_file)
            ewrite(0,*) "   Events  : ", n_events
            if (event_period_s > 0.0) then
                ewrite(0,*) "   Period  : ", event_period_s, " s"
            else
                ewrite(0,*) "   Mode    : non-periodic (absolute times)"
            end if
            ewrite(0,*) "=================================================="
            ewrite(0,*)

            ! Warn if the user is running with fixed dt : their stated dt will be overridden at every event. Better they know upfront
            if (.not. have_option('/timestepping/adaptive_timestep') .and. &
                .not. have_option('/solver_options/Non_Linear_Solver/Fixed_Point_Iteration/adaptive_timestep_nonlinear')) then
                ewrite(0,*)
                ewrite(0,*) "WARNING simulation_events: No adaptive timestepping detected. Your fixed timestep will be reduced whenever an event is approached. Consider enabling adaptive timestepping so that dt can recover automatically after each event OR set your fixed dt to be a divisor of all event intervals."
                ewrite(0,*)
            end if
        end if

    end subroutine load_simulation_events


    !> @brief Apply event-driven timestep constraint.
    !>
    !> Called at the start of each timestep. If the proposed dt would overshoot the
    !> next event, dt is reduced to land exactly on it.
    !> A warning is printed to the terminal ONLY if the event can't be reached
    !> because dt_min prevents snapping to it
    !>
    !> @param current_time  Simulation time at the START of this timestep (seconds).
    !> @param dt            Proposed timestep (seconds). Modified in place if an event
    !>                      falls within the current window
    !> @param dt_min        Optional minimum allowable timestep (seconds). A warning is
    !>                      printed if snapping would require going below this value
    subroutine apply_event_timestep(current_time, dt, dt_min)
        implicit none

        real, intent(in)           :: current_time
        real, intent(inout)        :: dt
        real, intent(in), optional :: dt_min

        integer :: i
        real    :: t_local, t_event_local, time_to_event, effective_dt_min

        if (.not. events_active) return
        if (n_events == 0)       return

        ! Apply post-event dt cap from previous timestep
        if (just_landed_on_event) then
            if (post_event_dt > 0.0) dt = min(dt, post_event_dt)
            just_landed_on_event = .false.
        end if

        effective_dt_min = 0.0
        if (present(dt_min)) effective_dt_min = dt_min

        ! Map current time into [0, period) for periodic schedules
        if (event_period_s > 0.0) then
            t_local = mod(current_time, event_period_s)
        else
            t_local = current_time
        end if

        do i = 1, n_events

            if (event_period_s > 0.0) then
                t_event_local = mod(event_times_s(i), event_period_s)
                time_to_event = t_event_local - t_local
                ! wrap negative differences around the period
                if (time_to_event < -EVENT_TOL) time_to_event = time_to_event + event_period_s
            else
                t_event_local = event_times_s(i)
                time_to_event = t_event_local - t_local
                ! For non periodic, skip events already in the past
                if (time_to_event < -EVENT_TOL) cycle
            end if

            ! Already sitting on this event ... nothing to do
            if (abs(time_to_event) <= EVENT_TOL) cycle

            ! Event falls within the current dt window
            if (time_to_event < dt - EVENT_TOL) then
                if (time_to_event < effective_dt_min - EVENT_TOL) then
                    ! Cannot snap: dt_min prevents reaching this event ... warn and skip
                    if (getprocno() == 1) then
                        ewrite(0,*) "WARNING simulation_events: cannot land on event at t =", &
                            current_time + time_to_event, "s"
                        ewrite(0,*) "  time_to_event =", time_to_event, &
                            "s is smaller than dt_min =", effective_dt_min, "s"
                        ewrite(0,*) "  Reduce dt_min to fix this."
                    end if
                else
                    dt = min(dt, time_to_event)
                    just_landed_on_event = .true.
                end if
                return
            end if

        end do

    end subroutine apply_event_timestep

end module multi_events