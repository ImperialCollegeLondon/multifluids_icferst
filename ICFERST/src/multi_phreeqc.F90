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


    subroutine init_PHREEQC(Mdims, packed_state, id, concetration_phreeqc)
        implicit none
        integer, INTENT(out) :: id
        integer :: i,j,k
        integer :: nxyz
        integer :: nthreads
        integer :: status
        integer :: save_on
  
        double precision, dimension(:), allocatable, target :: hydraulic_K
        double precision, dimension(:), allocatable   :: rv
        double precision, dimension(:), allocatable   :: por
        double precision, dimension(:), allocatable   :: sat
        character(100)                                :: string
        integer                                       :: ncomps, nspecies
        character(100),   dimension(:), allocatable   :: components, species_name
        integer,          dimension(:,:), allocatable :: ic1, ic2
        double precision, dimension(:,:), allocatable :: f1
        double precision, dimension(:,:), allocatable, INTENT(OUT) :: concetration_phreeqc
        double precision                              :: time, time_step
        double precision                              :: pH
        double precision, dimension(:,:), allocatable :: species_c
        type( state_type ), intent( inout ) :: packed_state
        type(multi_dimensions), intent( in ) :: Mdims
        integer :: cv_nod, iphase
        type(tensor_field), pointer :: H_field,O_field, charge_field, Ca_field, Cl_field
        type(tensor_field), pointer :: K_field, N_field, Na_field
  
        H_field=>extract_tensor_field(packed_state, "PackedPassiveTracer_H")
        O_field=>extract_tensor_field(packed_state, "PackedPassiveTracer_O")
        charge_field=>extract_tensor_field(packed_state, "PackedPassiveTracer_charge")
        Ca_field=>extract_tensor_field(packed_state, "PackedPassiveTracer_Ca")
        Cl_field=>extract_tensor_field(packed_state, "PackedPassiveTracer_Cl")
        K_field=>extract_tensor_field(packed_state, "PackedPassiveTracer_K")
        N_field=>extract_tensor_field(packed_state, "PackedPassiveTracer_N")
        Na_field=>extract_tensor_field(packed_state, "PackedPassiveTracer_Na")
  
        nxyz = Mdims%cv_nonods
        nthreads = 0
        id = RM_Create(nxyz, nthreads)
        status = RM_LoadDatabase(id, "phreeqc.dat.in")
  
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
        allocate(rv(nxyz))
        rv = 1.0
        status = RM_SetRepresentativeVolume(id, rv)
        ! Set initial porosity
        allocate(por(nxyz))
        por = 0.2
        status = RM_SetPorosity(id, por)
  
        !Read input file
        status = RM_RunFile(id, 1, 1, 1, "test.pqi.in")
        ! Clear contents of workers and utility
        string = "DELETE; -all"
        status = RM_RunString(id, 1, 0, 1, string)  ! workers, initial_phreeqc, utility
        ! Determine number of components to transport
        ncomps = RM_FindComponents(id)
  
        allocate(components(ncomps))
        do i = 1, ncomps
          status = RM_GetComponent(id, i, components(i))
        enddo
        !Setting up intial conditions - I think this is basically saying what kind
        !of reactions are expecting Phreeqc to run
        allocate(ic1(nxyz,7), ic2(nxyz,7), f1(nxyz,7))
        ic1 = -1
        ic2 = -1
        f1 = 1.0
        do i = 1, nxyz
          ic1(i,1) = 1       ! Solution 1
          ic1(i,2) = -1      ! Equilibrium phases 1
          ic1(i,3) = 1       ! Exchange none
          ic1(i,4) = -1      ! Surface none
          ic1(i,5) = -1      ! Gas phase none
          ic1(i,6) = -1      ! Solid solutions none
          ic1(i,7) = -1      ! Kinetics none
        enddo
  
        !Transfer solutions and reactants from the InitialPhreeqc instance
        !to the reaction-module workers - basically the part will actually be running things
        status = RM_InitialPhreeqc2Module(id, ic1, ic2, f1)
        time = 0.0
        time_step = 720.
        allocate(concetration_phreeqc(nxyz, ncomps))
        status = RM_SetTime(id, time)
        status = RM_SetTimeStep(id, time_step)
        save_on = RM_SetSpeciesSaveOn(id, 1)
  
        status = RM_RunCells(id)
        !Get the output data
        status = RM_GetConcentrations(id, concetration_phreeqc)
        H_field%val(1,1,:) = concetration_phreeqc(:,1)
        O_field%val(1,1,:) = concetration_phreeqc(:,2)
        charge_field%val(1,1,:) = concetration_phreeqc(:,3)
        Ca_field%val(1,1,:) = concetration_phreeqc(:,4)
        Cl_field%val(1,1,:) = concetration_phreeqc(:,5)
        K_field%val(1,1,:) = concetration_phreeqc(:,6)
        N_field%val(1,1,:) = concetration_phreeqc(:,7)
        Na_field%val(1,1,:) = concetration_phreeqc(:,8)
  
  
    !   print *, maxval(K_field%val)
  
      end subroutine init_PHREEQC
  
      subroutine testing_PHREEQC(Mdims, packed_state, id, concetration_phreeqc)
        implicit none
        integer , INTENT(INOUT) :: id
        integer :: status
  
        double precision, dimension(:,:),INTENT(INOUT) :: concetration_phreeqc
        double precision :: time, time_step
        type( state_type ), intent( inout ) :: packed_state
        type(multi_dimensions), intent( in ) :: Mdims
        type(tensor_field), pointer :: H_field,O_field, charge_field, Ca_field, Cl_field
        type(tensor_field), pointer :: K_field, N_field, Na_field
        integer :: cv_nod, iphase
  
        H_field=>extract_tensor_field(packed_state, "PackedPassiveTracer_H")
        O_field=>extract_tensor_field(packed_state, "PackedPassiveTracer_O")
        charge_field=>extract_tensor_field(packed_state, "PackedPassiveTracer_charge")
        Ca_field=>extract_tensor_field(packed_state, "PackedPassiveTracer_Ca")
        Cl_field=>extract_tensor_field(packed_state, "PackedPassiveTracer_Cl")
        K_field=>extract_tensor_field(packed_state, "PackedPassiveTracer_K")
        N_field=>extract_tensor_field(packed_state, "PackedPassiveTracer_N")
        Na_field=>extract_tensor_field(packed_state, "PackedPassiveTracer_Na")
  
  
        concetration_phreeqc(:,1) =  H_field%val(1,1,:)
        concetration_phreeqc(:,2) =  O_field%val(1,1,:)
        concetration_phreeqc(:,3) =  charge_field%val(1,1,:)
        concetration_phreeqc(:,4) =  Ca_field%val(1,1,:)
        concetration_phreeqc(:,5) =  Cl_field%val(1,1,:)
        concetration_phreeqc(:,6) =  K_field%val(1,1,:)
        concetration_phreeqc(:,7) =  N_field%val(1,1,:)
        concetration_phreeqc(:,8) =  Na_field%val(1,1,:)
        status = RM_SetConcentrations(id, concetration_phreeqc)
        status = RM_RunCells(id)
        !Get the output data
        status = RM_GetConcentrations(id, concetration_phreeqc)
        H_field%val(1,1,:) = concetration_phreeqc(:,1)
        O_field%val(1,1,:) = concetration_phreeqc(:,2)
        charge_field%val(1,1,:) = concetration_phreeqc(:,3)
        Ca_field%val(1,1,:) = concetration_phreeqc(:,4)
        Cl_field%val(1,1,:) = concetration_phreeqc(:,5)
        K_field%val(1,1,:) = concetration_phreeqc(:,6)
        N_field%val(1,1,:) = concetration_phreeqc(:,7)
        Na_field%val(1,1,:) = concetration_phreeqc(:,8)
  
  
      end subroutine testing_PHREEQC

end module multi_phreeqc