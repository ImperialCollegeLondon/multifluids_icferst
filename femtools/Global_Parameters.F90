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

module global_parameters
  !!< This routine exists to save us all from argument list hell!
  !!<
  !!< All the global parameters which don't change while fluidity is running
  !!< should live here. I am building this up as I encounter more parameters
  !!< in the code. It would be great if others did the same.
  !!<
  !!< The correct syntax for accessing this module is:
  !!<
  !!< use global_parameters, only: parameter1, parameter2 ...
  !!<
  !!< Try to only use the parameters which are needed locally.

  ! Debug specific paramaters are contained in fldebug_parameters
  ! (to resolve build dependencies)
  use fldebug_parameters
  use iso_c_binding

  implicit none

  !------------------------------------------------------------------------
  ! Precision parameters
  !------------------------------------------------------------------------
  !! Number of digits past the decimal point for a real
  integer, parameter :: real_digits_10 = precision(0.0)

  !! Integer size in bytes
  integer, parameter :: integer_size=bit_size(0)/8
  ! The real_size depends on reals having the obvious ieee754 sizes.
  !! Real size in bytes
#ifdef DOUBLEP
  integer, parameter :: real_size=8
#else
  integer, parameter :: real_size=4
#endif

  !------------------------------------------------------------------------
  ! Parameters controlling the scheme used in the flow core.
  !------------------------------------------------------------------------


  !! The simulation start time (model time)
  real, save :: simulation_start_time
  !! The simulation start CPU time
  real, save :: simulation_start_cpu_time
  !! The simulation start wall time
  real, save :: simulation_start_wall_time

  !! Accumulated system time.
  real, save, target :: current_time
  !! The timestep.
  real, save, target ::  dt
  !! The current timestep number
  integer, save, target :: timestep = 0

  real, parameter:: pi = 3.1415926535897931

  !------------------------------------------------------------------------
  ! Parameters for parallel
  !------------------------------------------------------------------------

  !! When upscaling a problem (e.g. from 16 to 32 processors),
  !! we want to run sam on 32 processors to do the domain decomposition.
  !! But only 16 processors will have data on disk. However,
  !! all 32 processors still have to go through populate_state
  !! to make sure it goes through all the MPI calls and doesn't
  !! deadlock. So we record whether the process is an "active" process,
  !! one that has data on disk.
  logical :: is_active_process = .true.
  integer :: no_active_processes = -1

  !------------------------------------------------------------------------
  ! Field names and paths
  !------------------------------------------------------------------------

  ! Zeroing long strings is EXPENSIVE.
  ! (See the commit message for r11059)
  ! That is why we supply an empty_name and an empty_path
  ! as, e.g.,
  ! field%option_path=empty_path
  ! is much much quicker than
  ! field%option_path="" .
  ! This is probably a bug in gcc, but it is a bug in gcc
  ! that we have to live with.

  !! Field names are permitted to be as long as Fortran names.
  integer, parameter :: FIELD_NAME_LEN=101
  character(len=FIELD_NAME_LEN) :: empty_name=""

  !! Maximum length of an option path
  !! -asiri - added an IFDEF statment for Intel Compiler, where a max char length is 7198 !!-ao
#ifndef using_gfortran
  integer, parameter :: OPTION_PATH_LEN=8192
  character(len=OPTION_PATH_LEN) :: empty_path=""

  !! Maximum length of a python string representing a function
  integer, parameter :: PYTHON_FUNC_LEN=8192
  character(len=PYTHON_FUNC_LEN) :: empty_python_func=""
#else
  integer, parameter :: OPTION_PATH_LEN=7198
  character(len=OPTION_PATH_LEN) :: empty_path=""

  !! Maximum length of a python string representing a function
  integer, parameter :: PYTHON_FUNC_LEN=7198
  character(len=PYTHON_FUNC_LEN) :: empty_python_func=""
#endif


  !! Name of the topology mesh in state - this mesh is used by adaptivity
  !! for the error metric etc.
  character(len=FIELD_NAME_LEN):: topology_mesh_name=""

  !! Name of mesh to be handled by adapt_state()
  character(len=FIELD_NAME_LEN):: adaptivity_mesh_name=""

  !! optionpath where the periodic boundary conditions are defined
  character(len=OPTION_PATH_LEN), dimension(3) :: periodic_boundary_option_path=""

  !! The bounding box of the input domain
  ! dim x 2
  ! (:, 1) are the minima along each coordinate
  ! (:, 2) are the maxima along each coordinate
  real, dimension(:, :), allocatable :: domain_bbox

  real :: domain_volume

  !! Are we using a new mesh?
  logical :: new_mesh = .false.
  !! New timestep for limiter
  logical :: new_lim = .false.
  !! When on-the-sphere, the planet radius is needed.
  ! The variable is initiliased as unity, to avoid garbage
  ! being passed around.
  real :: surface_radius = 1.0

  ! Colouring "enum".  These can't be in the colouring module due to circular dependencies
  integer, parameter :: COLOURING_CG1 = 1
  integer, parameter :: COLOURING_DG0 = 2
  integer, parameter :: COLOURING_DG1 = 3
  integer, parameter :: COLOURING_DG2 = 4
  integer, parameter :: NUM_COLOURINGS = 4

  !! Multiphase prototype, reservoir simulator
  logical :: is_porous_media = .false.
  !! Porous media model initialisation (Gravity capillary equilibration)
  logical :: is_porous_initialisation = .false.
  !! Multiphase prototype, models temperature
  logical :: has_temperature = .false.
  !! Mutliphase prototype, solve Stokes equations for momentum
  logical :: solve_stokes
  !! Multiphase prototype, models concentration
  logical :: has_concentration = .false.
  !! Checking if using the P0DG for velocity,
  !! special because we need to avoid the use of PressureMesh_Continuous
  logical :: is_P0DGP1 = .false.
  !! Has boussinesq option on
  logical :: has_boussinesq_aprox = .false.
  !! Porous media, anisotropic permeability
  logical :: has_anisotropic_permeability = .false.
  
  !!Public variable to be used in Adaptive_NonLinear to re-scale the effective convergence
  real :: backtrack_or_convergence
  logical :: FPI_have_converged = .false.
  logical :: after_adapt = .false.
  logical :: first_time_step = .false.
  logical :: first_nonlinear_time_step = .false.
  logical :: solver_not_converged = .false.
  !!Public string containing a generic warning and tips to get the code working
  character(len=OPTION_PATH_LEN) :: multi_generic_warning =""

contains

  function get_surface_radius() bind(c)
    !C-inter-operable subroutine for making the value of surface_radius availabe
    ! to C functions.
    implicit none

    real(kind=c_double) :: get_surface_radius

    get_surface_radius = real(surface_radius, kind=c_double)
  end function get_surface_radius

end module global_parameters
