!    Copyright (C) 2020 Imperial College London and others.
!
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Affero General Public License
!    as published by the Free Software Foundation,
!    version 3.0 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without seven the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA
#include "fdebug.h"

module multi_phreeqc

#ifdef USING_PHREEQC
  use IPhreeqc
  use PhreeqcRM
#endif

    use fldebug
    use state_module
    use fields
    use field_options
    use spud
    use parallel_tools
    use global_parameters
    use futils, only: int2str
    use boundary_conditions_from_options
    use parallel_tools, only : allmax, allmin, isparallel, getprocno
    use parallel_fields
    use memory_diagnostics
    use initialise_fields_module, only: initialise_field_over_regions
    use halos
    use multi_data_types
    use multi_tools

    implicit none

  contains

    subroutine init_PHREEQC(Mdims, packed_state, id, concetration_phreeqc, after_adapt)
        implicit none

        integer, INTENT(out) :: id
        type( state_type ), intent( inout ) :: packed_state
        type(multi_dimensions), intent( in ) :: Mdims
        real, dimension(:,:), allocatable, INTENT(OUT) :: concetration_phreeqc
        logical, intent (in) :: after_adapt

        integer :: i,j,k, ICncomp, iphase, icomp, cv_inod
        integer :: nxyz
        integer :: nthreads
        integer :: status
        integer :: save_on

        real, dimension(:), allocatable, target :: hydraulic_K
        real, dimension(:), allocatable   :: rv
        real, dimension(:), allocatable   :: sat
        character(100)                                :: string
        integer                                       :: ncomps, nspecies
        character(100),   dimension(:), allocatable   :: components, species_name
        integer,          dimension(:,:), allocatable :: ic1, ic2
        real, dimension(:,:), allocatable :: f1
        real  :: time, time_step
        real  :: pH
        real, dimension(:,:), allocatable :: species_c
        type(tensor_field), pointer :: tfield, tfield_old
        type(vector_field), pointer :: vfield
        character( len = option_path_len ) :: option_path, option_name
        character(len=option_path_len), dimension(:),  allocatable :: file_strings
        character(len=25), dimension(7) :: reaction_types
#ifdef USING_PHREEQC
        nxyz = Mdims%cv_nonods
        nthreads = 0
        id = RM_Create(nxyz, nthreads)
        call get_option("/porous_media/Phreeqc_coupling/database_name",option_name,default='phreeqc.dat.in')
        status = RM_LoadDatabase(id, trim(option_name))

        ! Set properties
        status = RM_SetErrorOn(id, 1)
        status = RM_SetErrorHandlerMode(id, 2)  ! exit on error
        status = RM_SetComponentH2O(id, 0)
        status = RM_SetRebalanceFraction(id, 0.5d0)
        status = RM_SetRebalanceByCell(id, 1)
        status = RM_UseSolutionDensityVolume(id, 0)
        status = RM_SetPartitionUZSolids(id, 0)

        ! Set concentration units
        status = RM_SetUnitsSolution(id, 2)      ! 1, mg/L; 2, mol/L; 3, kg/kgs
        status = RM_SetUnitsPPassemblage(id, 1)  ! 0, mol/L cell; 1, mol/L water; 2 mol/L rock
        status = RM_SetUnitsExchange(id, 1)      ! 0, mol/L cell; 1, mol/L water; 2 mol/L rock
        status = RM_SetUnitsSurface(id, 1)       ! 0, mol/L cell; 1, mol/L water; 2 mol/L rock
        status = RM_SetUnitsGasPhase(id, 1)      ! 0, mol/L cell; 1, mol/L water; 2 mol/L rock
        status = RM_SetUnitsSSassemblage(id, 1)  ! 0, mol/L cell; 1, mol/L water; 2 mol/L rock
        status = RM_SetUnitsKinetics(id, 1)      ! 0, mol/L cell; 1, mol/L water; 2 mol/L rock

        ! Set conversion from seconds to user units (days)
        !status = RM_SetTimeConversion(id, dble(1.0 / 86400.0))

        ! Set representative volume
        allocate(rv(nxyz)); rv = 1.0
        status = RM_SetRepresentativeVolume(id, rv)
        deallocate(rv)
        ! Set initial porosity
        vfield=>extract_vector_field(packed_state,"MeanPoreCV")
        status = RM_SetPorosity(id, vfield%val(1,:))

        !Read input file
        call get_option("/porous_media/Phreeqc_coupling/simulation_name",option_name)
        status = RM_RunFile(id, 1, 1, 1, trim(option_name))
        ! Clear contents of workers and utility
        string = "DELETE; -all"
        status = RM_RunString(id, 1, 0, 1, string)  ! workers, initial_phreeqc, utility
        ! Determine number of components to transport
        ncomps = RM_FindComponents(id)
        !Check how many species have been defined for Fotran and check that they match!
        ICncomp= 0
        do k = 1, option_count("/material_phase[0]/scalar_field")!We check the first phase only
          call get_option("/material_phase[0]/scalar_field["// int2str( k - 1)//"]/name", option_name)
          if (option_name(1:7)=="Species") then
            ICncomp = ICncomp + 1
          end if
        end do
        !PHREEQC will have an extra one which is charge
        if (ncomps /= ICncomp) then
          FLAbort("The number of Species defined in ICFERST and PHREEQC do not match. Please double check both input files.")
        end if

        allocate(components(ncomps))
        do i = 1, ncomps
          status = RM_GetComponent(id, i, components(i))
        enddo
        !Setting up intial conditions - I think this is basically saying what kind
        !of reactions are expecting Phreeqc to run
        reaction_types(1) = 'SOLUTION'
        reaction_types(2) = 'EQUILIBRIUM_PHASES'
        reaction_types(3) = 'EXCHANGE'
        reaction_types(4) = 'SURFACE'
        reaction_types(5) = 'GAS_PHASE'
        reaction_types(6) = 'SOLID_SOLUTIONS'
        reaction_types(7) = 'KINETICS'

        call read_inputfile(file_strings)

        allocate(ic1(nxyz,7), ic2(nxyz,7), f1(nxyz,7))
        ic1 = -1
        ic2 = -1
        f1 = 1.

       do k =1,1000
         if (len(trim(file_strings(k))) == 0) cycle
         do i = 1,7
           if (index(file_strings(k), trim(reaction_types(i))) /= 0) then
             do j = 1,nxyz
               ic1(j,i) = 1
             end do
           end if
         end do
       end do

        !Transfer solutions and reactants from the InitialPhreeqc instance
        !to the reaction-module workers - basically the part will actually be running things
        status = RM_InitialPhreeqc2Module(id, ic1, ic2, f1)
        call get_option( '/timestepping/timestep', time_step )
        call get_option( '/timestepping/current_time', time )
        allocate(concetration_phreeqc(nxyz, ncomps))
        status = RM_SetTime(id, time)
        status = RM_SetTimeStep(id, time_step)
        save_on = RM_SetSpeciesSaveOn(id, 1)
        if (time <time_step+1.) then
          status = RM_RunCells(id)
          if (.not. after_adapt) then
          !Get the output data
          status = RM_GetConcentrations(id, concetration_phreeqc)
            do iphase = 1, 1!Mdims%nphase!TODO SINGLE PHASE FOR THE TIME BEING, WE NEED TO KNOW HOW PHREEQC WOULD DEAL WITH MULTIPHASE
              do icomp = 1, ncomps
                tfield=>extract_tensor_field(packed_state,get_packed_Species_name(components(icomp), .false.))
                tfield_old=>extract_tensor_field(packed_state,get_packed_Species_name(components(icomp), .true.))
                do cv_inod = 1, Mdims%cv_nonods!Since PHREEQC is not following column major, we use a do loop to hopefully speed it up
                  tfield%val(1,iphase,cv_inod) = concetration_phreeqc(cv_inod, icomp)
                  tfield_old%val(1,iphase,cv_inod) = concetration_phreeqc(cv_inod, icomp)
                end do
              end do
            end do
          end if
        end if
#endif
      end subroutine init_PHREEQC

      subroutine run_PHREEQC(Mdims, packed_state, id, concetration_phreeqc)
        implicit none
        integer , INTENT(INOUT) :: id
        integer :: status, iphase, icomp, cv_inod, ncomps, i

        real, dimension(:,:),INTENT(INOUT) :: concetration_phreeqc
        real :: time, time_step
        type( state_type ), intent( inout ) :: packed_state
        type(multi_dimensions), intent( in ) :: Mdims
        type(tensor_field), pointer :: tfield
#ifdef USING_PHREEQC        
        ! Determine number of components to transport
        ncomps = RM_FindComponents(id)
        allocate(components(ncomps))
        do i = 1, ncomps
          status = RM_GetComponent(id, i, components(i))
        enddo

        do iphase = 1, 1!Mdims%nphase!SINGLE PHASE FOR THE TIME BEING, WE NEED TO KNOW HOW PHREEQC WOULD DEAL WITH MULTIPHASE
          do icomp = 1, ncomps
            tfield=>extract_tensor_field(packed_state,get_packed_Species_name(components(icomp), .false.))
            do cv_inod = 1, Mdims%cv_nonods!Since PHREEQC is not following column major, we use a do loop to hopefully speed it up
              concetration_phreeqc(cv_inod, icomp) = tfield%val(1,iphase,cv_inod)
            end do
          end do
        end do

        status = RM_SetConcentrations(id, concetration_phreeqc)
        status = RM_RunCells(id)
        !Get the output data
        status = RM_GetConcentrations(id, concetration_phreeqc)
        do iphase = 1, 1!Mdims%nphase!SINGLE PHASE FOR THE TIME BEING, WE NEED TO KNOW HOW PHREEQC WOULD DEAL WITH MULTIPHASE
          do icomp = 1, ncomps
            tfield=>extract_tensor_field(packed_state,get_packed_Species_name(components(icomp), .false.))
            do cv_inod = 1, Mdims%cv_nonods!Since PHREEQC is not following column major, we use a do loop to hopefully speed it up
              tfield%val(1,iphase,cv_inod) = concetration_phreeqc(cv_inod, icomp)
            end do
          end do
        end do
      end subroutine run_PHREEQC

      subroutine read_inputfile(file_strings)

        implicit none

        integer :: i, start, ierr
        character( len = option_path_len ) :: option_path, option_name
        character(len=option_path_len), dimension(:),  allocatable, intent(out) :: file_strings

        allocate(file_strings(1000))
        call get_option("/porous_media/Phreeqc_coupling/simulation_name",option_name)
        open(unit= 89, file=trim(option_name), status='old', action='read')
        i = 1
        start = 1
        do while (.true.)
            read(89,'(a)', iostat=ierr) file_strings(i)
            if (ierr/=0) exit
            i = i + 1
        end do
        close(89)

#endif
      end subroutine

      subroutine deallocate_PHREEQC(id)

        implicit none

        integer, intent(in) :: id
        integer :: status

        status = RM_Destroy(id)
        
      end subroutine


  !>@author Geraldine Regnier, Pablo Salinas
  !>@brief: Finds the field name in diamond given a name in PHREEQC. Fields have the convention of being named in ICFERST as SPECIES_component,
  !> for example Species_O for oxygen.
    function get_packed_Species_name(PHREEQC_name, old_field)
      implicit none
      character(len = *), INTENT(IN) :: PHREEQC_name
      logical, intent(in) :: old_field
      character (len = option_path_len) :: get_packed_Species_name
      !Local variables
      integer :: k, buffer
      character( len = option_path_len ) :: option_name

      do k = 1, option_count("/material_phase[0]/scalar_field")!We check the first phase only
        call get_option("/material_phase[0]/scalar_field["// int2str( k - 1)//"]/name", option_name)
        if (option_name(1:7)=="Species") then
          !To avoid finding names that are not supposed to be, for example if we have Carbon and Calcite,
          !we would have C and Ca we dont want to identify Ca as Carbon we limit the lenght of the search
          !Check that the name contains the PHREEQC name, if it does, then use that name
          buffer = 8+len(trim(PHREEQC_name))
          if (option_name(8:) == "_"//trim(PHREEQC_name)) then
            if (old_field) then
              get_packed_Species_name = "PackedOldSpecies"//trim(option_name(8:))
            else
              get_packed_Species_name = "PackedSpecies"//trim(option_name(8:))
            end if
          end if
        end if
      end do
    end function

end module multi_phreeqc
