  
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
  

module multi_data_types
    use fldebug
    use sparse_tools_petsc
    use sparse_tools
    use fields_data_types
    use fields_allocates
    use global_parameters
    use state_module
    use fields
    use spud
    use multi_tools


    interface allocate_multi_dev_shape_funs
        module procedure allocate_multi_dev_shape_funs1
        module procedure allocate_multi_dev_shape_funs2
        module procedure allocate_multi_dev_shape_funs3
    end interface

    interface allocate_multi_field
        module procedure allocate_multi_field1
        module procedure allocate_multi_field2
    end interface

    type multi_dimensions
        integer :: ndim       !Number of dimensions
        integer :: cv_nloc    !Number of local control volumes
        integer :: u_nloc     !Number of local velocity nodes
        integer :: cv_snloc   !Number of local control volumes on the surface?
        integer :: u_snloc    !Number of local velocity nodes on the surface?
        integer :: nstate     !Number of states in state
        integer :: ncomp      !Number of components
        integer :: xu_nloc    !Number of local velocity nodes of the Continuous mesh
        integer :: x_nloc     !Number of local control volumes of the Continuous mesh
        integer :: x_snloc    !Number of local surface control volumes of the Continuous mesh
        integer :: x_nloc_p1  !???
        integer :: x_nonods_p1!???
        integer :: p_nloc     !Number of local pressure nodes
        integer :: p_snloc    !Number of local pressure nodes on the surface?
        integer :: mat_nloc   !Number of local material nodes
        integer :: totele     !Total number of elements
        integer :: stotel     !Total number of surface elements?
        integer :: cv_nonods  !Total number of control volumes
        integer :: p_nonods   !Total number of pressure nodes
        integer :: mat_nonods !Total number of sub-control volumes
        integer :: u_nonods   !Total number of velocity nodes
        integer :: xu_nonods  !Total number of velocity nodes of the Continuous mesh
        integer :: x_nonods   !Total number of control volumes of the Continuous mesh
        integer :: ph_nloc    !Number of ????
        integer :: ph_nonods  !Total number of ????
        integer :: nphase     !Total number of phases
        integer :: npres      !Total number of pressure
        integer :: n_in_pres  !nphase/npres

    end type multi_dimensions


    type multi_discretization_opts
        !This type includes the necessary information to choose from the different discretization options available
        integer ::  cv_ele_type, p_ele_type, u_ele_type, mat_ele_type, u_sele_type, cv_sele_type
        integer ::  t_disopt, v_disopt
        real ::     t_beta, v_beta, t_theta, v_theta, u_theta
        integer ::  t_dg_vel_int_opt, u_dg_vel_int_opt, v_dg_vel_int_opt, w_dg_vel_int_opt, &
                    in_ele_upwind, dg_ele_upwind, nits_flux_lim_t, nits_flux_lim_volfra, nits_flux_lim_comp
        logical ::  volfra_use_theta_flux, volfra_get_theta_flux, comp_use_theta_flux, &
                    comp_get_theta_flux, t_use_theta_flux, t_get_theta_flux, scale_momentum_by_volume_fraction
    end type multi_discretization_opts

    type multi_gi_dimensions
        integer :: cv_ngi     !Number of gauss integer points
        integer :: scvngi     !Number of gauss integer points in the surface of a control volume
        integer :: sbcvngi    !Number of gauss integer points in the surface boundary of a control volume
        integer :: nface      !Number of faces per element
    end type multi_gi_dimensions


    !Data structure to store all the shape functions to facilitate its movement throughtout the code
    type multi_shape_funs
        real, pointer, dimension( : , : ) :: cvn => null()! Control volume shape function; dimension( cv_nloc, cv_ngi )
        real, pointer, dimension( : ) :: cvweight=> null()!Weigth of the control volume; dimension( cv_ngi )
        real, pointer, dimension(  : , : ) :: cvfen=> null()!Finite element of the control volume; dimension( cv_nloc, cv_ngi )
        real, pointer, dimension( : , : , :)  ::  cvfenlx_all=> null()!dimension( ndim, cv_nloc, cv_ngi )
        real, pointer, dimension(  : , : )  :: ufen=> null()!Finite element of the element; dimension( u_nloc, cv_ngi )
        real, pointer, dimension(  : , :,: )  :: ufenlx_all=> null()!dimension( ndim, u_nloc, cv_ngi )
        integer, pointer, dimension(  : , : )  :: cv_neiloc=> null()!dimension( cv_nloc, scvngi )
        logical, pointer, dimension(  : , : ) :: cv_on_face=> null(), cvfem_on_face=> null()!dimension( cv_nloc, scvngi )
        real, pointer, dimension(  : , : )  :: scvfen=> null(), scvfenslx=> null(), scvfensly=> null()!dimension( cv_nloc, scvngi )
        real, pointer, dimension( : )  :: scvfeweigh=> null()!dimension( scvngi )
        real, pointer, dimension(  : , : ,: )  :: scvfenlx_all=> null()!dimension( ndim, cv_nloc, scvngi )
        real, pointer, dimension(  : , : )  :: sufen=> null(), sufenslx=> null(), sufensly=> null()!dimension( u_nloc, scvngi )
        real, pointer, dimension(  : , :, : )  :: sufenlx_all=> null() !dimension( ndim, u_nloc, scvngi )
        logical, pointer, dimension(  : , : )  :: u_on_face=> null(), ufem_on_face=> null()!dimension( u_nloc, scvngi )
        real, pointer, dimension(  : , : )  :: sbcvn=> null()!dimension( cv_snloc, sbcvngi )
        real, pointer, dimension(  : , : )  :: sbcvfen=> null(), sbcvfenslx=> null(), sbcvfensly=> null()! dimension( cv_snloc, sbcvngi )
        real, pointer, dimension( : )  :: sbcvfeweigh=> null()!dimension( sbcvngi )
        real, pointer, dimension(  : , :, : )  :: sbcvfenlx_all=> null()!dimension( ndim, cv_snloc, sbcvngi )
        real, pointer, dimension(  : , : )  :: sbufen=> null(), sbufenslx=> null(), sbufensly=> null()!dimension( u_snloc, sbcvngi )
        real, pointer, dimension(  : , :,: )  :: sbufenlx_all=> null() !dimension( ndim, u_snloc, sbcvngi )
        integer, pointer, dimension(  : , : )  :: cv_sloclist=> null()!dimension( nface, cv_snloc )
        integer, pointer, dimension(  : , : )  :: u_sloclist=> null()!dimension( nface, u_snloc )
        integer, pointer, dimension( : )  :: findgpts=> null()!dimension( cv_nloc + 1 )
        integer, pointer, dimension( : )  :: colgpts=> null()!dimension( cv_nloc * scvngi )
        integer :: ncolgpts
        type(petsc_csr_matrix) ::CV2FE !Matrix to convert from CV to FE
        type(petsc_csr_matrix) ::FE2CV !Matrix to convert from FE to CV

    end type multi_shape_funs

    !Data structure to store the derivatives of the shape functions and conversors from reference element to local
    type multi_dev_shape_funs
        real :: volume!Volume of the local element
        real, pointer, dimension(:) :: detwei=> null()!Determinant times weigth (i.e: conversor from reference element to local element)
        real, pointer, dimension(:) :: ra => null()   !???
        real, pointer, dimension(:, :, :) :: cvfenx_all=> null()!Space derivatives of the pressure (CV) shape functions
        real, pointer, dimension(:, :, :) :: ufenx_all=> null()!Space derivatives of the velocity (FE) shape functions
        real, pointer, dimension(:, :, :) :: nx_all => null()!Space derivatives of a generic field.
        real, pointer, dimension(:, :, :) :: inv_jac => null()!Inverse of the Jacobian matrix
    end type multi_dev_shape_funs
    !This type comprises the four necessary variables to represent matrices using a CSR structure
    type multi_sparsity
        integer :: ncol
        integer, pointer, dimension(:) :: fin=> null()
        integer, pointer, dimension(:) :: col=> null()
        integer, pointer, dimension(:) :: mid=> null()
    end type multi_sparsity

    !This data type contains all the sparsities necessary in the multiphase prototype code
    type multi_sparsities
        type (multi_sparsity) :: acv     !CV multi-phase eqns (e.g. vol frac, temp)
        type (multi_sparsity) :: small_acv !Local CV multi-phase eqns (e.g. vol frac, temp)
        type (multi_sparsity) :: ele     !Element connectivity
        type (multi_sparsity) :: dgm_pha !Force balance sparsity
        type (multi_sparsity) :: ct      !CT sparsity - global continuity eqn
        type (multi_sparsity) :: C       !C sparsity operating on pressure in force balance
        type (multi_sparsity) :: cmc     !pressure matrix for projection method
        type (multi_sparsity) :: m       !CV-FEM matrix
        type (multi_sparsity) :: ph      !ph matrix
    end type multi_sparsities

    type multi_ndgln
        integer, dimension( : ), pointer  :: cv=> null()     !Control volume local to global numbering
        integer, dimension( : ), pointer  :: u=> null()      !Velocity local to global numbering
        integer, dimension( : ), pointer  :: p=> null()      !Pressure local to global numbering
        integer, dimension( : ), pointer  :: x=> null()      !Continuous mesh pressure local to global numbering
        integer, dimension( : ), pointer  :: x_p1=> null()   !Continuous mesh pressure P1 local to global numbering
        integer, dimension( : ), pointer  :: xu=> null()     !Continuous mesh velocity local to global numbering
        integer, dimension( : ), pointer  :: mat=> null()    !Pressure discontinuous local to global numbering
        integer, dimension( : ), pointer ::  suf_cv=> null() !Surface control volume local to global numering
        integer, dimension( : ), pointer ::  suf_p=> null()  !Pressure local to global numbering
        integer, dimension( : ), pointer ::  suf_u=> null()  !Velocity surface local to global numbering
    end type multi_ndgln

    type multi_matrices
        real, dimension( :, :, : ), pointer :: C => null()!Pressure matrix using a FE discretization (storable)
        real, dimension( :, :, : ), pointer :: C_CV => null()!Pressure matrix using a CV discretization (storable)
        REAL, DIMENSION( :, :, : ), pointer :: U_RHS => null()!Rigth hand side of the momentum equation
        real, dimension( :, :, : ), pointer :: CT => null()!Continuity equation matrix
        type(vector_field) :: CT_RHS!Rigth hand side of the continuity equation
        type(petsc_csr_matrix) :: petsc_ACV!Matrix of the saturation equation
        type(vector_field) :: CV_RHS!Rigth hand side of the saturation equation
        real, dimension( :, :, : ), pointer :: PIVIT_MAT => null()!Mass matrix (matrix form by the sigmas) (storable)
        integer, dimension(:), pointer :: ICOLOR => null()!Array used to accelerate the creation of CMC in COLOR_GET_CMC_PHA_FAST
        integer :: NCOLOR !Number of colors in ICOLOR
        type(petsc_csr_matrix):: DGM_PETSC!Big matrix to solve the pressure in inertia flows (don't know much more)
        logical :: NO_MATRIX_STORE !Flag to whether calculate and use DGM_PETSC or C
        logical :: CV_pressure     !Flag to whether calculate the pressure using FE (ASSEMB_FORCE_CTY) or CV (cv_assemb)
        logical :: stored = .false.!Flag to be true when the storable matrices have been stored
        logical :: compact_PIVIT_MAT = .false. !Flag to know whether to use a compacted mass matrix or not
    end type multi_matrices


    type porous_adv_coefs
        real, dimension( :, :, :, : ), pointer :: adv_coef => null()!Sigmas at the boundary to calculate fluxes
        real, dimension( :, :, :, : ), pointer :: inv_adv_coef => null()!Inverse of sigmas at the boundary to calculate fluxes
        real, dimension( :, :, :, : ), pointer :: adv_coef_grad => null()!Gradient of the sigmas at the boundary to calculate fluxes
    end type porous_adv_coefs

    type multi_field
        real, dimension( :, :, :, : ), pointer :: val => null()

        logical :: have_field = .false. ! do we need this field for this simulation?
        logical :: is_constant = .false. ! if ( .true. ) nonods = 1 for what follows     -   DELETE THIS MAYBE ???

        integer :: memory_type = -1 ! 0 Isotropic tensor - ( 1, 1, nphase, nonods ) - this is unrolled as ( ndim, ndim, nphase, nonods )
                                    ! 1 Isotropic - ( 1, 1, nphase, nonods ) - diagonal
                                    ! 2 Anisotropic - ( ndim, ndim, nphase, nonods )
                                    ! 3 Isotropic coupled - ( 1, nphase, nphase, nonods )
                                    ! 4 Anisotropic coupled (aka Full Metal Jacket) - ( 1, ndim x nphase, ndim x nphase, nonods )

        integer :: ndim1 = -1, ndim2 = -1, ndim3 = -1 ! dimensions of field

    end type multi_field


    type multi_absorption
        !Comprises all the absorption terms that migth be required
        type (multi_field) :: PorousMedia ! <= Always memory_type = 2
        type (multi_field) :: Components
        type (multi_field) :: Temperature
        type (multi_field) :: Velocity
        type (multi_field) :: Flooding
    end type multi_absorption

    type multi_transport_scalar
        !Contains all the information for generic scalar fields to solve for
        character(len = FIELD_NAME_LEN) :: name!To extract from state
        character(len = option_path_len) :: path!Path from diamond
        logical :: coupled_field!Is the field coupled between phases?
        type (multi_field) :: absorption!Absorption of this field
    end type

    type multi_pipe_package
        !Contains all the information required to model pipes.
        real, dimension( :, :, : ), pointer  :: gamma_pres_abs=> null()
        real, dimension( :, :, : ), pointer  :: gamma_pres_abs_nano=> null()
        real, dimension( : ), pointer        :: mass_pipe=> null()
        real, dimension( : ), pointer        :: mass_cvfem2pipe=> null()
        real, dimension( : ), pointer        :: mass_pipe2cvfem=> null()
        real, dimension( : ), pointer        :: mass_cvfem2pipe_true=> null()
    end type

    type multi_outfluxes
        !Contains variables to analyse the flux across the BCs that the user is interested
        logical :: calculate_flux !True if all the process related with this has to start or not
        integer, dimension(:), allocatable :: outlet_id !ids the user wants
        real :: porevolume ! for outfluxes.csv to calculate the pore volume injected
        real, allocatable, dimension(:,:,:) :: totout!(field -saturation, temperature-, Mdims%nphase, size(outlet_id))
        real, dimension(:,:),  allocatable  :: intflux
    end type

     type pipe_coords
            integer :: ele, npipes                               !Element containing pipes, pipes per element
            logical, allocatable, dimension(:) :: pipe_index     !nodes with pipes
            integer, allocatable, dimension(:) :: pipe_corner_nds1!size npipes
            integer, allocatable, dimension(:) :: pipe_corner_nds2!size npipes
     end type pipe_coords

    private :: allocate_multi_dev_shape_funs1, allocate_multi_dev_shape_funs2, allocate_multi_dev_shape_funs3,&
         allocate_multi_field1, allocate_multi_field2

contains

    subroutine allocate_multi_field1( state, Mdims, field_name, mfield )
        !*********UNTESTED*********
        implicit none

        type( state_type ), intent( in ) :: state
        type( multi_dimensions ), intent(in) :: Mdims
        character( len = * ), intent( in ) :: field_name

        type( multi_field ), intent( inout ) :: mfield

        type( tensor_field ), pointer :: tfield
        integer :: ndim, nphase, nonods, stat, dimensions

        mfield%have_field = .true.

        tfield => extract_tensor_field( state, trim( field_name ), stat )

        ndim = Mdims%ndim ; nphase = Mdims%nphase

        if ( stat /= 0 ) FLAbort( "Cannot determine multi_field source." )

        ! Decide whether the field is constant throught the domain or not
        if ( have_option(trim(tfield%option_path)//"prescribed") ) mfield%is_constant = .true.  ! This logic is not correct

        ! Number of nodes of the field
        nonods = size( tfield%val, 3 )

        ! Number of dimensions of the coupling, for example ndim*ndim*nphase
        call get_option( trim(tfield%option_path) // "/type/dimensions", dimensions, default = -1)
        if ( dimensions < 1 ) FLAbort( "Wrong input for dimensions." )

        ! Select memory type
        if ( have_option(trim(tfield%option_path) // "/type/Anisotropic_coupled" ) ) then
            mfield%memory_type = 4
        else if ( have_option(trim(tfield%option_path) // "/type/Anisotropic" ) ) then
            mfield%memory_type = 2
        else if ( have_option(trim(tfield%option_path) // "/type/Isotropic_coupled" ) ) then
            mfield%memory_type = 3
        else if ( have_option(trim(tfield%option_path) // "/type/Isotropic" ) ) then
            mfield%memory_type = 1
            if ( trim( tfield%name ) == "Viscosity" ) mfield%memory_type = 0
        else
            FLAbort( "Wrong memory type selected." )
        end if

        select case ( mfield%memory_type )
            case( 0, 1 ) ! Isotropic ( full and diagonal )
                mfield%ndim1 = 1    ; mfield%ndim2 = 1           ; mfield%ndim3 = nphase
            case( 2 )    ! Anisotropic
                mfield%ndim1 = ndim ; mfield%ndim2 = ndim        ; mfield%ndim3 = nphase
            case( 3 )    ! Isotropic coupled
                mfield%ndim1 = 1    ; mfield%ndim2 = nphase      ; mfield%ndim3 = nphase
            case( 4 )    ! Anisotropic coupled
                mfield%ndim1 = 1    ; mfield%ndim2 = ndim*nphase ; mfield%ndim3 = ndim*nphase
            case default
                FLAbort( "Cannot determine multi_field memrory_type." )
        end select

        mfield%val(1:mfield%ndim1, 1:mfield%ndim2, 1:mfield%ndim3, 1:nonods) => tfield%val

        return
    end subroutine allocate_multi_field1

    subroutine allocate_multi_field2( Mdims, mfield, nonods_in, field_name)
        !*********UNTESTED*********
        !PorousMedia_AbsorptionTerm tested
        implicit none
        integer, intent(in) :: nonods_in!Number of nodes of the field.
        type( multi_dimensions ), intent(in) :: Mdims
        type( multi_field ), intent( inout ) :: mfield
        character( len = * ), optional, intent( in ) :: field_name
        !Local variables
        integer :: ndim, nphase, nonods

        mfield%have_field = .true.

        ndim = Mdims%ndim ; nphase = Mdims%nphase

        !Number of nodes of the field
        nonods = nonods_in
        !Decide whether the field is constant throught the domain or not
        mfield%is_constant = (nonods == 1)

        if (present(field_name)) then
            !Depending on the field, different possibilities
            if (trim(field_name)=="PorousMedia_AbsorptionTerm") then
                mfield%is_constant = .false.!For porous media it cannot be constant
                mfield%memory_type = 2 !We force this memory despite not being the most comprised
                !because it enables us to remove copies of memory and because for real 3D problems it is very unlikely that
                !the permeability will be isotropic in all the regions
            end if

            if (trim(field_name)=="Flooding_AbsorptionTerm") then
                mfield%is_constant = .false.!It cannot be constant
                mfield%memory_type = 1!The absorption is always isotropic
            end if

            if (trim(field_name)=="ComponentAbsorption") then
                mfield%memory_type = 4
                nonods = Mdims%cv_nonods
                mfield%is_constant = .false.
            end if
            if (trim(field_name)=="SourceTerm") then
                mfield%memory_type = 5
                nonods = nonods_in
                mfield%is_constant = .false.
            end if
        end if
        select case ( mfield%memory_type )
            case( 0, 1 ) ! Isotropic ( full and diagonal )
                mfield%ndim1 = 1    ; mfield%ndim2 = 1           ; mfield%ndim3 = nphase
            case( 2 )    ! Anisotropic
                mfield%ndim1 = ndim ; mfield%ndim2 = ndim        ; mfield%ndim3 = nphase
            case( 3 )    ! Isotropic coupled
                mfield%ndim1 = 1    ; mfield%ndim2 = nphase      ; mfield%ndim3 = nphase
            case( 4 )    ! Anisotropic coupled
                mfield%ndim1 = 1    ; mfield%ndim2 = ndim*nphase ; mfield%ndim3 = ndim*nphase
            case ( 5 )   !Source term
                mfield%ndim1 = ndim    ; mfield%ndim2 = nphase ; mfield%ndim3 = 1
            case default
                FLAbort( "Cannot determine multi_field memrory_type." )
        end select
        !Allocate and initialize memory
        allocate(mfield%val(1:mfield%ndim1, 1:mfield%ndim2, 1:mfield%ndim3, 1:nonods)); mfield%val = 0.

        return
    end subroutine allocate_multi_field2



    subroutine deallocate_multi_field(mfield, and_destroy)
        implicit none
        type( multi_field ), intent( inout ) :: mfield
        logical, optional, intent(in) :: and_destroy

        if(present_and_true(and_destroy)) then
            deallocate(mfield%val)
        end if

        nullify(mfield%val)
        mfield%memory_type = 0
        mfield%have_field = .false.
        mfield%is_constant = .false.
        mfield%memory_type = -1
        mfield%ndim1 = -1; mfield%ndim2 = -1; mfield%ndim3 = -1

    end subroutine deallocate_multi_field

    subroutine deallocate_multi_absorption(multi_absorp, and_destroy)
        implicit none
        type( multi_absorption ), intent( inout ) :: multi_absorp
        logical, optional, intent(in) :: and_destroy
        !Local variables
        logical :: and_destroy2

        and_destroy2 = present_and_true(and_destroy)

        if (associated(multi_absorp%PorousMedia%val)) call deallocate_multi_field(multi_absorp%PorousMedia, and_destroy2)
        if (associated(multi_absorp%Components%val))  call deallocate_multi_field(multi_absorp%Components, and_destroy2)
        if (associated(multi_absorp%Temperature%val)) call deallocate_multi_field(multi_absorp%Temperature, and_destroy2)
        if (associated(multi_absorp%Velocity%val))    call deallocate_multi_field(multi_absorp%Velocity, and_destroy2)
        if (associated(multi_absorp%Flooding%val))    call deallocate_multi_field(multi_absorp%Flooding, and_destroy2)

    end subroutine deallocate_multi_absorption

    subroutine get_multi_field(mfield, inode_in, output)
        implicit none
        integer, intent(in) :: inode_in
        type( multi_field ), intent( in ) :: mfield
        real, dimension(:,:),intent( inout ) :: output!must have size(ndim*nphase, ndim*nphase)
        !local variables
        integer :: iphase, jphase, idim, jdim, ndim, inode

        inode = inode_in
        if (mfield%is_constant) inode = 1
        select case (mfield%memory_type)
            case (0)!Isotropic full
                ndim = size(output,2)/mfield%ndim3
                do iphase = 1, mfield%ndim3!nphase
                    output(1 + (iphase-1)*ndim:ndim + (iphase-1)*ndim, 1 + (iphase-1)*ndim:ndim + (iphase-1)*ndim) =&
                         mfield%val(1,1,iphase,inode)
                end do
            case (1)!Isotropic
                output = 0.;ndim = size(output,2)/mfield%ndim3
                do iphase = 1, mfield%ndim3!nphase
                    do idim = 1, ndim
                        output(idim+(iphase-1)*ndim, idim+(iphase-1)*ndim) = mfield%val(1,1,iphase,inode)
                    end do
                end do
            case (2)!Anisotropic
                output = 0.
                do iphase = 1, mfield%ndim3!nphase
                    do jdim = 1, mfield%ndim2!ndim
                        do idim = 1, mfield%ndim2!ndim
                            output(idim+(iphase-1)*mfield%ndim2,jdim+(iphase-1)*mfield%ndim2) = mfield%val(idim,jdim,iphase,inode)
                        end do
                    end do
                end do
            case (3)!isotropic coupled
                output = 0.;ndim = size(output,2)/mfield%ndim3
                do iphase = 1, mfield%ndim3
                    do jphase = 1, mfield%ndim3
                        do idim = 1, ndim
                            output(idim+(iphase-1)*ndim ,idim+(jphase-1)*ndim) = mfield%val(1,iphase,jphase,inode)
                        end do
                    end do
                end do
            case default !Anisotropic coupled
                output = mfield%val(1,:,:,inode)
        end select

    end subroutine get_multi_field


    subroutine get_multi_field_inverse(mfield, inode_in, output)
        implicit none
        integer, intent(in) :: inode_in
        type( multi_field ), intent( in ) :: mfield
        real, dimension(:,:),intent( inout ) :: output!must have size(ndim*nphase, ndim*nphase)
        !local variables
        integer :: iphase, jphase, idim, ndim, inode

        inode = inode_in
        if (mfield%is_constant) inode = 1
        select case (mfield%memory_type)
            case (1)!Isotropic
                output = 0.;ndim = size(output,2)/mfield%ndim3
                do iphase = 1, mfield%ndim3!nphase
                    do idim = 1, ndim
                        output(idim+(iphase-1)*ndim, idim+(iphase-1)*ndim) = 1./mfield%val(1,1,iphase,inode)
                    end do
                end do
            case (2)!Anisotropic
                output = 0.
                do iphase = 1, mfield%ndim3!nphase
                    output(1+(iphase-1)*mfield%ndim2:mfield%ndim2+(iphase-1)*mfield%ndim2,&
                        1+(iphase-1)*mfield%ndim2:mfield%ndim2+(iphase-1)*mfield%ndim2) = &
                                inverse(mfield%val(:,:,iphase,inode))
                end do
            case (3)!isotropic coupled
                output = 0.;ndim = size(output,2)/mfield%ndim3
                do iphase = 1, mfield%ndim3
                    do jphase = 1, mfield%ndim3
                        do idim = 1, ndim
                            output(idim+(iphase-1)*ndim, idim+(jphase-1)*ndim) = 1./mfield%val(1,iphase,jphase,inode)
                        end do
                    end do
                end do
            case default !Anisotropic coupled
                output = inverse(mfield%val(1,:,:,inode))
        end select



    end subroutine get_multi_field_inverse

    subroutine print_multi_field(mfield, inode_in, dimension)
        implicit none
        integer, intent(in) :: inode_in
        type( multi_field ), intent( in ) :: mfield
        integer :: dimension!must have size(ndim*nphase, ndim*nphase)
        !Local variables
        real, dimension(dimension,dimension) :: Matrix

        call get_multi_field(mfield, inode_in, Matrix)
        call printMatrix(Matrix)

    end subroutine print_multi_field

    subroutine add_array_to_multi_field(mfield, b, xpos, ypos, inode)
        !mfield = mfield + b
        !xpos and ypos are the starting positions
        !for a full matrix they have to be one
        implicit none
        integer, intent(in) :: inode
        integer, intent(in) :: xpos, ypos
        type( multi_field ), intent( inout ) :: mfield
        real, dimension(:,:), intent(in) :: b
        !Local variables
        integer :: fxpos, fypos, idim, jdim, iphase, ndim, jphase

        fxpos = xpos + size(b,1) - 1
        fypos = ypos + size(b,2) - 1
        ndim = size(b,2)/mfield%ndim3
        select case (mfield%memory_type)
            case (0)!Isotropic viscosity
                do iphase = 1 + (xpos-1)/ndim, fxpos/ndim!1, mfield%ndim3
                    mfield%val(1,1,iphase,inode) = mfield%val(1,1,iphase,inode) + &
                    b(1+(xpos-1)/ndim+(iphase-1)*ndim,1+(ypos-1)/ndim+(iphase-1)*ndim)
                end do
            case (1)!Isotropic
                do iphase = 1 + (xpos-1)/ndim, fxpos/ndim!1, mfield%ndim3
                    mfield%val(1,1,iphase,inode) = mfield%val(1,1,iphase,inode) + &
                        b(1+(iphase-1)*ndim ,1+(iphase-1)*ndim)
                end do
            case (2)!Anisotropic
                !Work out the involved phases from the position
                do iphase = 1 + (xpos-1)/ndim, fxpos/ndim!!do iphase =1,mfield%ndim3
                    do idim = 1 + (xpos-1)/iphase, fxpos/mfield%ndim3!jdim = 1, mfield%ndim2!ndim
                        do jdim = 1 + (xpos-1)/iphase, fypos/mfield%ndim3!idim = 1, mfield%ndim2!ndim
                            mfield%val(idim,jdim,iphase,inode) = mfield%val(idim,jdim,iphase,inode) +&
                                 b(idim+(iphase-1)*mfield%ndim2,jdim+(iphase-1)*mfield%ndim2)
                        end do
                    end do
                end do
            case (3)!isotropic coupled
                do iphase = 1 + (xpos-1)/ndim, fxpos/ndim!1, mfield%ndim3
                    do jphase = 1 + (ypos-1)/ndim, fypos/ndim!1, mfield%ndim3
                        mfield%val(1,iphase,jphase,inode) = mfield%val(1,iphase,jphase,inode) +&
                         b(1+(iphase-1)*ndim ,1+(jphase-1)*ndim)
                    end do
                end do
            case default !Anisotropic coupled
                mfield%val(1,xpos:fxpos,ypos:fypos,inode) = mfield%val(1,xpos:fxpos,ypos:fypos,inode) + b(xpos:fxpos,ypos:fypos)
        end select

    end subroutine add_array_to_multi_field

    subroutine add_multi_field_to_array(mfield, b, xpos, ypos, inode, a_in)
        !b = b + a * mfield
        !xpos and ypos are the starting positions
        !for a full matrix they have to be one
        implicit none
        integer, intent(in) :: inode
        integer, intent(in) :: xpos, ypos
        type( multi_field ), intent( in ) :: mfield
        real, dimension(:,:), intent(inout) :: b
        real, optional, intent(in) :: a_in
        !Local variables
        integer :: fxpos, fypos, idim, jdim, iphase, ndim, jphase
        real :: a

        if(present(a_in)) then
            a = a_in
        else
            a = 1.0
        end if
        fxpos = xpos + size(b,1) - 1
        fypos = ypos + size(b,2) - 1
        ndim = size(b,2)/mfield%ndim3
        select case (mfield%memory_type)
            case (0)!Isotropic viscosity
                do iphase = 1 + (xpos-1)/ndim, fxpos/ndim!1, mfield%ndim3
                    b(1+(xpos-1)/ndim+(iphase-1)*ndim: fxpos/mfield%ndim3 +(iphase-1)*ndim,1+(ypos-1)/ndim+(iphase-1)*ndim:fxpos/mfield%ndim3 +(iphase-1)*ndim)=&
                    b(1+(xpos-1)/ndim+(iphase-1)*ndim: fxpos/mfield%ndim3 +(iphase-1)*ndim,1+(ypos-1)/ndim+(iphase-1)*ndim:fxpos/mfield%ndim3 +(iphase-1)*ndim)+&
                        + a * mfield%val(1,1,iphase,inode)
                end do
            case (1)!Isotropic
                do iphase = 1 + (xpos-1)/ndim, fxpos/ndim!1, mfield%ndim3
                    do idim = 1 + (xpos-1)/iphase, fxpos/mfield%ndim3
                        b(idim+(iphase-1)*ndim ,idim+(iphase-1)*ndim) = &
                            b(idim+(iphase-1)*ndim ,idim+(iphase-1)*ndim) + a * mfield%val(1,1,iphase,inode)
                    end do
                end do
            case (2)!Anisotropic
                !Work out the involved phases from the position
                do iphase = 1 + (xpos-1)/ndim, fxpos/ndim!!do iphase =1,mfield%ndim3
                    do idim = 1 + (xpos-1)/iphase, fxpos/mfield%ndim3!jdim = 1, mfield%ndim2!ndim
                        do jdim = 1 + (ypos-1)/iphase, fypos/mfield%ndim3!idim = 1, mfield%ndim2!ndim
                            b(idim+(iphase-1)*mfield%ndim2,jdim+(iphase-1)*mfield%ndim2) = &
                                b(idim+(iphase-1)*mfield%ndim2,jdim+(iphase-1)*mfield%ndim2) +a * mfield%val(idim,jdim,iphase,inode)
                        end do
                    end do
                end do
            case (3)!isotropic coupled
                do iphase = 1 + (xpos-1)/ndim, fxpos/ndim!1, mfield%ndim3
                    do jphase = 1 + (ypos-1)/ndim, fypos/ndim!1, mfield%ndim3
                        do idim = 1 + (xpos-1)/mfield%ndim3, fxpos/mfield%ndim3
                            b(idim+(iphase-1)*ndim ,idim+(jphase-1)*ndim) = &
                                b(idim+(iphase-1)*ndim ,idim+(jphase-1)*ndim) + a * mfield%val(1,iphase,jphase,inode)
                        end do
                    end do
                end do
            case default !Anisotropic coupled
                b(xpos:fxpos,ypos:fypos) = b(xpos:fxpos,ypos:fypos) + a * mfield%val(1,xpos:fxpos,ypos:fypos,inode)
        end select

    end subroutine add_multi_field_to_array

    subroutine mult_multi_field_by_array(mfield, b, inode)
        !mfield = mfield * b
        implicit none
        integer, intent(in) :: inode
        type( multi_field ), intent( inout ) :: mfield
        real, dimension(:,:), intent(in) :: b
        !Local variables
        integer ::  iphase, ndim, jphase
        real, dimension(:,:), allocatable :: miniB

        ndim = size(b,2)/mfield%ndim3
        select case (mfield%memory_type)
            case (0, 1)!Isotropic
                do iphase = 1, mfield%ndim3
                    mfield%val(1,1,iphase,inode) = mfield%val(1,1,iphase,inode) * &
                        b(1+(iphase-1)*ndim ,1+(iphase-1)*ndim)
                end do
            case (2)!Anisotropic
                !Work out the involved phases from the position
                do iphase =1,mfield%ndim3
                    mfield%val(:,:,iphase,inode) = matmul(mfield%val(:,:,iphase,inode), &
                            b(1+(iphase-1)*mfield%ndim2:iphase*mfield%ndim2 ,1+(iphase-1)*mfield%ndim2:iphase*mfield%ndim2))
                end do
            case (3)!isotropic coupled
                allocate(miniB(mfield%ndim3,mfield%ndim3))!(nphase,nphase)
                do iphase = 1, mfield%ndim3!nphase
                    do jphase = 1, mfield%ndim3!nphase
                        miniB(iphase, jphase) = b(1+(iphase-1)*ndim ,1+(jphase-1)*ndim)
                    end do
                end do
                mfield%val(1,:,:,inode) = matmul(mfield%val(1,:,:,inode),miniB)
                deallocate(miniB)
            case default !Anisotropic coupled
                mfield%val(1,:,:,inode) = matmul(mfield%val(1,:,:,inode),b)
        end select

    end subroutine mult_multi_field_by_array


    subroutine mult_multi_field_by_array_on_array(mfield, b, inode)
        !b = mfield * b
        implicit none
        integer, intent(in) :: inode
        type( multi_field ), intent( in ) :: mfield
        real, dimension(:,:), intent(inout) :: b
        !Local variables
        integer :: idim, jdim, iphase, ndim, jphase
        real, dimension(:,:), allocatable :: miniB

        ndim = size(b,2)/mfield%ndim3
        select case (mfield%memory_type)
            case (0)!Isotropic viscosity
                do iphase = 1, mfield%ndim3
                    b(1+(iphase-1)*ndim: ndim+(iphase-1)*ndim,1+(iphase-1)*ndim:ndim+(iphase-1)*ndim) = &
                        b(1+(iphase-1)*ndim: ndim+(iphase-1)*ndim,1+(iphase-1)*ndim:ndim+(iphase-1)*ndim) *&
                             mfield%val(1,1,iphase,inode)
                end do
            case (1)!Isotropic
                do iphase = 1, mfield%ndim3
                    do idim = 1, ndim
                        b(idim+(iphase-1)*ndim ,idim+(iphase-1)*ndim) = &
                            b(idim+(iphase-1)*ndim ,idim+(iphase-1)*ndim) * mfield%val(1,1,iphase,inode)
                    end do
                end do
            case (2)!Anisotropic
                ndim = size(b,2)/mfield%ndim3
                !Work out the involved phases from the position
                do iphase =1,mfield%ndim3
                    b(1+(iphase-1)*mfield%ndim2:iphase*mfield%ndim2 ,1+(iphase-1)*mfield%ndim2:iphase*mfield%ndim2) = &
                    matmul(mfield%val(:,:,iphase,inode), b(1+(iphase-1)*mfield%ndim2:iphase*mfield%ndim2 ,1+(iphase-1)*mfield%ndim2:iphase*mfield%ndim2))
                end do
            case (3)!Isotropic coupled
                ndim = size(b,2)/mfield%ndim3
                allocate(miniB(mfield%ndim3,mfield%ndim3))!(nphase,nphase)
                do iphase = 1, mfield%ndim3!nphase
                    do jphase = 1, mfield%ndim3!nphase
                        miniB(iphase, jphase) = b(1+(iphase-1)*ndim ,1+(jphase-1)*ndim)
                    end do
                end do
                miniB = matmul(mfield%val(1,:,:,inode),miniB)
                do iphase = 1, mfield%ndim3!nphase
                    do jphase = 1, mfield%ndim3!nphase
                        do idim = 1 , ndim
                            b(idim+(iphase-1)*ndim ,idim+(jphase-1)*ndim) = miniB(iphase, jphase)
                        end do
                    end do
                end do
                deallocate(miniB)
            case default !Anisotropic coupled
                b = matmul(mfield%val(1,:,:,inode),b)
        end select

    end subroutine mult_multi_field_by_array_on_array



    subroutine scale_multi_field(mfield, a, inode)
        !mfield = a * mfield
        implicit none
        integer, intent(in) :: inode
        type( multi_field ), intent( inout ) :: mfield
        real, intent(in) :: a

        mfield%val(:,:,:,inode) = a * mfield%val(:,:,:,inode)

    end subroutine scale_multi_field

    subroutine linearise_multi_field( mfield, Mdims, ndgln )
      !*********UNTESTED*********
      implicit none
      type( multi_field ), intent( inout ) :: mfield
      type( multi_dimensions ), intent( in ) :: Mdims
      integer, dimension( : ), pointer, intent( in ) :: ndgln

      integer, dimension( : ), pointer :: nodes
      integer :: ndim, nloc, ele

      ndim = Mdims%ndim ; nloc = Mdims%Mat_nloc ! This nloc should be improved in the future

      if ( nloc/=6 .or. nloc/=10 ) FLAbort( "I can only linearise P2 fields..." )

      do ele = 1, Mdims%totele
         nodes => ndgln( (ele-1)*nloc+1 : ele*nloc )

         mfield%val(:,:,:,nodes(2))=0.5*(mfield%val(:,:,:,nodes(1))+mfield%val(:,:,:,nodes(3)))
         mfield%val(:,:,:,nodes(4))=0.5*(mfield%val(:,:,:,nodes(1))+mfield%val(:,:,:,nodes(6)))
         mfield%val(:,:,:,nodes(5))=0.5*(mfield%val(:,:,:,nodes(3))+mfield%val(:,:,:,nodes(6)))

         if ( ndim > 2 ) then
            mfield%val(:,:,:,nodes(7))=0.5*(mfield%val(:,:,:,nodes(1))+mfield%val(:,:,:,nodes(10)))
            mfield%val(:,:,:,nodes(8))=0.5*(mfield%val(:,:,:,nodes(3))+mfield%val(:,:,:,nodes(10)))
            mfield%val(:,:,:,nodes(9))=0.5*(mfield%val(:,:,:,nodes(6))+mfield%val(:,:,:,nodes(10)))
         end if
      end do

      return
    end subroutine linearise_multi_field

















    subroutine allocate_multi_shape_funs(shape_fun,  Mdims, GIdims)
    !This subroutine allocates all the arrays in a multi_shape_funs data type
        implicit none
        type(multi_shape_funs), intent(inout) :: shape_fun
        type(multi_dimensions), intent(in) :: Mdims
        type(multi_GI_dimensions), intent(in) :: GIdims

        !Proceed to allocate the variables
        allocate(shape_fun%cvn(Mdims%cv_nloc, GIdims%cv_ngi) )
        allocate(shape_fun%cvweight(GIdims%cv_ngi))
        allocate(shape_fun%cvfen(Mdims%cv_nloc, GIdims%cv_ngi ))
        allocate(shape_fun%cvfenlx_all( Mdims%ndim, Mdims%cv_nloc, GIdims%cv_ngi ))
        allocate(shape_fun%ufen( Mdims%u_nloc, GIdims%cv_ngi ))
        allocate(shape_fun%ufenlx_all( Mdims%ndim, Mdims%u_nloc, GIdims%cv_ngi ))
        allocate(shape_fun%cv_neiloc( Mdims%cv_nloc, GIdims%scvngi ))
        allocate(shape_fun%cv_on_face( Mdims%cv_nloc, GIdims%scvngi ))
        allocate(shape_fun%cvfem_on_face( Mdims%cv_nloc, GIdims%scvngi ))
        allocate(shape_fun%scvfen( Mdims%cv_nloc, GIdims%scvngi ))
        allocate(shape_fun%scvfenslx( Mdims%cv_nloc, GIdims%scvngi ))
        allocate(shape_fun%scvfensly( Mdims%cv_nloc, GIdims%scvngi ))
        allocate(shape_fun%scvfeweigh( GIdims%scvngi ))
        allocate(shape_fun%scvfenlx_all(Mdims%ndim, Mdims%cv_nloc, GIdims%scvngi))
        allocate(shape_fun%sufen( Mdims%u_nloc, GIdims%scvngi ))
        allocate(shape_fun%sufenslx( Mdims%u_nloc, GIdims%scvngi ))
        allocate(shape_fun%sufensly( Mdims%u_nloc, GIdims%scvngi ))
        allocate(shape_fun%sufenlx_all( Mdims%ndim, Mdims%u_nloc, GIdims%scvngi ))
        allocate(shape_fun%u_on_face( Mdims%u_nloc, GIdims%scvngi ))
        allocate(shape_fun%ufem_on_face( Mdims%u_nloc, GIdims%scvngi ))
        allocate(shape_fun%sbcvn(Mdims%cv_snloc, GIdims%sbcvngi))
        allocate(shape_fun%sbcvfen( Mdims%cv_snloc, GIdims%sbcvngi ))
        allocate(shape_fun%sbcvfenslx( Mdims%cv_snloc, GIdims%sbcvngi ))
        allocate(shape_fun%sbcvfensly( Mdims%cv_snloc, GIdims%sbcvngi ))
        allocate(shape_fun%sbcvfeweigh(GIdims%sbcvngi))
        allocate(shape_fun%sbcvfenlx_all( Mdims%ndim, Mdims%cv_snloc, GIdims%sbcvngi ))
        allocate(shape_fun%sbufen(Mdims%u_snloc, GIdims%sbcvngi))
        allocate(shape_fun%sbufenslx(Mdims%u_snloc, GIdims%sbcvngi))
        allocate(shape_fun%sbufensly(Mdims%u_snloc, GIdims%sbcvngi))
        allocate(shape_fun%sbufenlx_all( Mdims%ndim, Mdims%u_snloc, GIdims%sbcvngi ))
        allocate(shape_fun%cv_sloclist( GIdims%nface, Mdims%cv_snloc ))
        allocate(shape_fun%u_sloclist( GIdims%nface, Mdims%u_snloc ))
        allocate(shape_fun%findgpts(Mdims%cv_nloc + 1))
        allocate(shape_fun%colgpts( Mdims%cv_nloc * GIdims%scvngi ))

        !Set everything to zero
        shape_fun%cvn = 0.
        shape_fun%cvweight = 0.
        shape_fun%cvfen = 0.
        shape_fun%cvfenlx_all = 0.
        shape_fun%ufen = 0.
        shape_fun%ufenlx_all = 0.
        shape_fun%cv_neiloc = 0
        shape_fun%cv_on_face = .false.
        shape_fun%cvfem_on_face = .false.
        shape_fun%scvfen = 0.
        shape_fun%scvfenslx = 0.
        shape_fun%scvfensly = 0.
        shape_fun%scvfeweigh = 0.
        shape_fun%scvfenlx_all = 0.
        shape_fun%sufen = 0.
        shape_fun%sufenslx = 0.
        shape_fun%sufensly = 0.
        shape_fun%sufenlx_all = 0.
        shape_fun%u_on_face = .false.
        shape_fun%ufem_on_face = .false.
        shape_fun%sbcvn = 0.
        shape_fun%sbcvfen = 0.
        shape_fun%sbcvfenslx = 0.
        shape_fun%sbcvfensly = 0.
        shape_fun%sbcvfeweigh = 0.
        shape_fun%sbcvfenlx_all = 0.
        shape_fun%sbufen = 0.
        shape_fun%sbufenslx = 0.
        shape_fun%sbufensly = 0.
        shape_fun%sbufenlx_all = 0.
        shape_fun%cv_sloclist = 0
        shape_fun%u_sloclist = 0
        shape_fun%findgpts = 0
        shape_fun%colgpts = 0

    end subroutine allocate_multi_shape_funs


    subroutine deallocate_multi_shape_funs(shape_fun)
        !This subroutine deallocates all the fields inside multi_shape_funs
        implicit none
        type(multi_shape_funs), intent(inout) :: shape_fun

        !Proceed to deallocate the variables
        deallocate(shape_fun%cvn)
        deallocate(shape_fun%cvweight)
        deallocate(shape_fun%cvfen)
        deallocate(shape_fun%cvfenlx_all)
        deallocate(shape_fun%ufen)
        deallocate(shape_fun%ufenlx_all)
        deallocate(shape_fun%cv_neiloc)
        deallocate(shape_fun%cv_on_face)
        deallocate(shape_fun%cvfem_on_face)
        deallocate(shape_fun%scvfen)
        deallocate(shape_fun%scvfenslx)
        deallocate(shape_fun%scvfensly)
        deallocate(shape_fun%scvfeweigh)
        deallocate(shape_fun%scvfenlx_all)
        deallocate(shape_fun%sufen)
        deallocate(shape_fun%sufenslx)
        deallocate(shape_fun%sufensly)
        deallocate(shape_fun%sufenlx_all)
        deallocate(shape_fun%u_on_face)
        deallocate(shape_fun%ufem_on_face)
        deallocate(shape_fun%sbcvn)
        deallocate(shape_fun%sbcvfen)
        deallocate(shape_fun%sbcvfenslx)
        deallocate(shape_fun%sbcvfensly)
        deallocate(shape_fun%sbcvfeweigh)
        deallocate(shape_fun%sbcvfenlx_all)
        deallocate(shape_fun%sbufen)
        deallocate(shape_fun%sbufenslx)
        deallocate(shape_fun%sbufensly)
        deallocate(shape_fun%sbufenlx_all)
        deallocate(shape_fun%cv_sloclist)
        deallocate(shape_fun%u_sloclist)
        deallocate(shape_fun%findgpts)
        deallocate(shape_fun%colgpts)
        if (associated(shape_fun%CV2FE%refcount)) call deallocate(shape_fun%CV2FE)
        if (associated(shape_fun%FE2CV%refcount)) call deallocate(shape_fun%FE2CV)

        !Nullify pointers
        nullify(shape_fun%cvn)
        nullify(shape_fun%cvweight)
        nullify(shape_fun%cvfen)
        nullify(shape_fun%cvfenlx_all)
        nullify(shape_fun%ufen)
        nullify(shape_fun%ufenlx_all)
        nullify(shape_fun%cv_neiloc)
        nullify(shape_fun%cv_on_face)
        nullify(shape_fun%cvfem_on_face)
        nullify(shape_fun%scvfen)
        nullify(shape_fun%scvfenslx)
        nullify(shape_fun%scvfensly)
        nullify(shape_fun%scvfeweigh)
        nullify(shape_fun%scvfenlx_all)
        nullify(shape_fun%sufen)
        nullify(shape_fun%sufenslx)
        nullify(shape_fun%sufensly)
        nullify(shape_fun%sufenlx_all)
        nullify(shape_fun%u_on_face)
        nullify(shape_fun%ufem_on_face)
        nullify(shape_fun%sbcvn)
        nullify(shape_fun%sbcvfen)
        nullify(shape_fun%sbcvfenslx)
        nullify(shape_fun%sbcvfensly)
        nullify(shape_fun%sbcvfeweigh)
        nullify(shape_fun%sbcvfenlx_all)
        nullify(shape_fun%sbufen)
        nullify(shape_fun%sbufenslx)
        nullify(shape_fun%sbufensly)
        nullify(shape_fun%sbufenlx_all)
        nullify(shape_fun%cv_sloclist)
        nullify(shape_fun%u_sloclist)
        nullify(shape_fun%findgpts)
        nullify(shape_fun%colgpts)
        nullify(shape_fun%FE2CV%refcount)
        nullify(shape_fun%CV2FE%refcount)

    end subroutine deallocate_multi_shape_funs

    subroutine deallocate_projection_matrices(shape_fun)
        !This subroutine deallocates projection matrices CV2FE & FE2CV stored in shape_fun
        implicit none
        type(multi_shape_funs), intent(inout) :: shape_fun

        if (associated(shape_fun%CV2FE%refcount)) call deallocate(shape_fun%CV2FE)
        if (associated(shape_fun%FE2CV%refcount)) call deallocate(shape_fun%FE2CV)
        nullify(shape_fun%FE2CV%refcount)
        nullify(shape_fun%CV2FE%refcount)

    end subroutine deallocate_projection_matrices

    !This subroutine, despite it can be called by itself it is highly recommended to be called ONLY through multi_sparsity/Get_Sparsity_Patterns
    subroutine allocate_multi_sparsities(Mspars, Mdims, mx_ncolacv, mx_ncoldgm_pha, mx_nct, mx_nc, mx_ncolm, mx_ncolph)
        !This subroutine allocates part of the memory inside Mspars
        implicit none
        type (multi_sparsities), intent(inout) :: Mspars
        type(multi_dimensions), intent(in) :: Mdims
        integer :: mx_ncolacv, mx_ncoldgm_pha, mx_nct, mx_nc, mx_ncolm, mx_ncolph

        if(.not.associated(Mspars%ACV%fin))          allocate( Mspars%ACV%fin( Mdims%cv_nonods * Mdims%nphase + 1 ))
        if(.not.associated(Mspars%ACV%col))          allocate( Mspars%ACV%col( mx_ncolacv ))
        if(.not.associated(Mspars%ACV%mid))          allocate(  Mspars%ACV%mid( Mdims%cv_nonods * Mdims%nphase ))

        if(.not.associated(Mspars%DGM_PHA%fin))      allocate(  Mspars%DGM_PHA%fin( Mdims%u_nonods * Mdims%nphase * Mdims%ndim + 1 ))
        if(.not.associated(Mspars%DGM_PHA%col))      allocate(  Mspars%DGM_PHA%col( mx_ncoldgm_pha ))
        if(.not.associated(Mspars%DGM_PHA%mid))      allocate(  Mspars%DGM_PHA%mid( Mdims%u_nonods * Mdims%nphase * Mdims%ndim ))

        if(.not.associated(Mspars%CT%fin))           allocate(  Mspars%CT%fin( Mdims%cv_nonods + 1 ))
        if(.not.associated(Mspars%CT%col))           allocate(  Mspars%CT%col( mx_nct ))
        if(.not.associated(Mspars%CT%mid))           allocate(  Mspars%CT%mid(  Mdims%cv_nonods ))

        if(.not.associated(Mspars%C%fin))            allocate(  Mspars%C%fin( Mdims%u_nonods + 1 ))
        if(.not.associated(Mspars%C%col))            allocate(  Mspars%C%col( mx_nc ))
        if(.not.associated(Mspars%C%mid))            allocate(  Mspars%C%mid( Mdims%u_nonods))

        if(.not.associated(Mspars%CMC%fin))          allocate( Mspars%CMC%fin( Mdims%cv_nonods + 1 ))
        if(.not.associated(Mspars%CMC%col))          allocate(  Mspars%CMC%col( 0 ))
        if(.not.associated(Mspars%CMC%mid))          allocate( Mspars%CMC%mid( Mdims%cv_nonods ))

        if(.not.associated( Mspars%M%fin))           allocate(  Mspars%M%fin( Mdims%cv_nonods + 1 ))
        if(.not.associated( Mspars%M%col))           allocate(  Mspars%M%col( mx_ncolm ))
        if(.not.associated( Mspars%M%mid))           allocate(  Mspars%M%mid( Mdims%cv_nonods ))

        if(.not.associated(Mspars%ph%fin))           allocate(  Mspars%ph%fin( Mdims%ph_nonods + 1 ))
        if(.not.associated(Mspars%ph%col))           allocate(  Mspars%ph%col( mx_ncolph ) )
        if(.not.associated(Mspars%ph%mid))           allocate(  Mspars%ph%mid( Mdims%ph_nonods ))

        Mspars%CT%col = 0 ; Mspars%C%fin = 0 ; Mspars%C%col = 0 ; Mspars%CMC%fin = 0
        Mspars%CMC%col = 0 ; Mspars%CMC%mid = 0 ; Mspars%M%fin = 0
        Mspars%M%col = 0 ; Mspars%M%mid = 0 ; Mspars%ph%fin = 0 ; Mspars%ph%col = 0

    end subroutine allocate_multi_sparsities

    subroutine deallocate_multi_sparsities(Mspars)
        !This subroutine deallocates all the memory inside Mspars
        implicit none
        type (multi_sparsities), intent(inout) :: Mspars

        !Proceed to deallocate sparsities
        if (associated(Mspars%acv%fin))       deallocate(Mspars%acv%fin)
        if (associated(Mspars%acv%col))       deallocate(Mspars%acv%col)
        if (associated(Mspars%acv%mid))       deallocate(Mspars%acv%mid)

        if (associated(Mspars%small_acv%fin)) deallocate(Mspars%small_acv%fin)
        if (associated(Mspars%small_acv%col)) deallocate(Mspars%small_acv%col)
        if (associated(Mspars%small_acv%mid)) deallocate(Mspars%small_acv%mid)

        if (associated(Mspars%dgm_pha%fin))   deallocate(Mspars%dgm_pha%fin)
        if (associated(Mspars%dgm_pha%col))   deallocate(Mspars%dgm_pha%col)
        if (associated(Mspars%dgm_pha%mid))   deallocate(Mspars%dgm_pha%mid)

        if (associated(Mspars%ct%fin))        deallocate(Mspars%ct%fin)
        if (associated(Mspars%ct%col))        deallocate(Mspars%ct%col)
        if (associated(Mspars%ct%mid))        deallocate(Mspars%ct%mid)

        if (associated(Mspars%C%fin))         deallocate(Mspars%C%fin)
        if (associated(Mspars%C%col))         deallocate(Mspars%C%col)
        if (associated(Mspars%C%mid))         deallocate(Mspars%C%mid)

        if (associated(Mspars%CMC%fin))       deallocate(Mspars%CMC%fin)
        if (associated(Mspars%CMC%col))       deallocate(Mspars%CMC%col)
        if (associated(Mspars%CMC%mid))       deallocate(Mspars%CMC%mid)

        if (associated(Mspars%M%fin))         deallocate(Mspars%M%fin)
        if (associated(Mspars%M%col))         deallocate(Mspars%M%col)
        if (associated(Mspars%M%mid))         deallocate(Mspars%M%mid)

        if (associated(Mspars%ph%fin))        deallocate(Mspars%ph%fin)
        if (associated(Mspars%ph%col))        deallocate(Mspars%ph%col)
        if (associated(Mspars%ph%mid))        deallocate(Mspars%ph%mid)

        !Proceed to nullify
        if (associated(Mspars%acv%fin))       nullify(Mspars%acv%fin)
        if (associated(Mspars%acv%col))       nullify(Mspars%acv%col)
        if (associated(Mspars%acv%mid))       nullify(Mspars%acv%mid)

        if (associated(Mspars%small_acv%fin)) nullify(Mspars%small_acv%fin)
        if (associated(Mspars%small_acv%col)) nullify(Mspars%small_acv%col)
        if (associated(Mspars%small_acv%mid)) nullify(Mspars%small_acv%mid)

        if (associated(Mspars%dgm_pha%fin))   nullify(Mspars%dgm_pha%fin)
        if (associated(Mspars%dgm_pha%col))   nullify(Mspars%dgm_pha%col)
        if (associated(Mspars%dgm_pha%mid))   nullify(Mspars%dgm_pha%mid)

        if (associated(Mspars%ct%fin))        nullify(Mspars%ct%fin)
        if (associated(Mspars%ct%col))        nullify(Mspars%ct%col)
        if (associated(Mspars%ct%mid))        nullify(Mspars%ct%mid)

        if (associated(Mspars%C%fin))         nullify(Mspars%C%fin)
        if (associated(Mspars%C%col))         nullify(Mspars%C%col)
        if (associated(Mspars%C%mid))         nullify(Mspars%C%mid)

        if (associated(Mspars%CMC%fin))       nullify(Mspars%CMC%fin)
        if (associated(Mspars%CMC%col))       nullify(Mspars%CMC%col)
        if (associated(Mspars%CMC%mid))       nullify(Mspars%CMC%mid)

        if (associated(Mspars%M%fin))         nullify(Mspars%M%fin)
        if (associated(Mspars%M%col))         nullify(Mspars%M%col)
        if (associated(Mspars%M%mid))         nullify(Mspars%M%mid)

        if (associated(Mspars%ph%fin))        nullify(Mspars%ph%fin)
        if (associated(Mspars%ph%col))        nullify(Mspars%ph%col)
        if (associated(Mspars%ph%mid))        nullify(Mspars%ph%mid)

    end subroutine deallocate_multi_sparsities

    subroutine allocate_multi_ndgln(ndgln, Mdims)
        implicit none
        type(multi_dimensions), intent(in) :: Mdims
        type(multi_ndgln), intent(inout) :: ndgln
        !Only allocate these three fields since the others are pointers to state
        allocate( ndgln%suf_cv( Mdims%stotel * Mdims%cv_snloc ), &
                  ndgln%suf_p( Mdims%stotel * Mdims%p_snloc ),&
                  ndgln%suf_u( Mdims%stotel * Mdims%u_snloc ) )
    end subroutine allocate_multi_ndgln

    subroutine deallocate_multi_ndgln(ndgln)
        implicit none
        type(multi_ndgln), intent(inout) :: ndgln
        !Only deallocate these three fields since the others are pointers to state
        deallocate( ndgln%suf_cv, ndgln%suf_p, ndgln%suf_u)
        nullify( ndgln%suf_cv, ndgln%suf_p, ndgln%suf_u)

    end subroutine deallocate_multi_ndgln

    subroutine destroy_multi_matrices(Mmat)
        !Deallocates all the memory used, intended to be used after mesh adapt
        implicit none
        type (multi_matrices), intent(inout) :: Mmat

        !Deallocate and nullify as required
        if (Mmat%CV_pressure) then!Deallocate just one of the two
            if (associated(Mmat%C_CV)) then
                deallocate (Mmat%C_CV); nullify(Mmat%C_CV)!this method won't be compatible with two matrices, but we want
            end if                                      !to get rid of it anyway (its unstable)
        else
            if (associated(Mmat%C)) then
                deallocate (Mmat%C); nullify(Mmat%C)
            end if
        end if
        if (associated(Mmat%U_RHS)) then
            deallocate (Mmat%U_RHS); nullify(Mmat%U_RHS)
        end if
        if (associated(Mmat%CT)) then
            deallocate (Mmat%CT); nullify(Mmat%CT)
        end if
        if (associated(Mmat%PIVIT_MAT)) then
            deallocate (Mmat%PIVIT_MAT); nullify(Mmat%PIVIT_MAT)
        end if
        if (associated(Mmat%CT_RHS%val)) call deallocate(Mmat%CT_RHS)
        !This one below gives problems unless it has been allocated at some point, as Mmat%CV_RHS%val is by default associated...
!        if (associated(Mmat%CV_RHS%val)) call deallocate(Mmat%CV_RHS)!<=Should not need to deallocate anyway as it is done somewhere else
!        if (associated(Mmat%petsc_ACV%refcount)) call deallocate(Mmat%petsc_ACV)!<=Should not need to deallocate anyway as it is done somewhere else
        if (associated(Mmat%DGM_PETSC%refcount)) call deallocate(Mmat%DGM_PETSC)
        !Set flag to recalculate
        Mmat%stored = .false.
    end subroutine destroy_multi_matrices

    subroutine allocate_multi_dev_shape_funs1(funs, DevFuns, nx_all_FE_size)
        implicit none
        type (multi_shape_funs), intent(in) ::funs
        type (multi_dev_shape_funs), intent(inout) :: DevFuns
        logical, optional, intent(in) :: nx_all_FE_size!If true then the size of the generic nx_all is set like
                                                       !like ufenlx_all, otherwise, like cvfenlx_all
        !Local variable
        logical :: nx_all_FE_size2
        nx_all_FE_size2 = .false.

        allocate(DevFuns%detwei(size(funs%cvfenlx_all,3)))
        allocate(DevFuns%ra(size(funs%cvfenlx_all,3)))
        allocate(DevFuns%cvfenx_all(size(funs%cvfenlx_all,1),size(funs%cvfenlx_all,2),size(funs%cvfenlx_all,3)))
        allocate(DevFuns%ufenx_all(size(funs%ufenlx_all,1),size(funs%ufenlx_all,2),size(funs%ufenlx_all,3)))
        allocate(DevFuns%inv_jac(size(funs%ufenlx_all,1),size(funs%ufenlx_all,1),size(funs%ufenlx_all,3)))

        if (present(nx_all_FE_size)) then
            if (nx_all_FE_size) nx_all_FE_size2 = .true.
        end if

        if (nx_all_FE_size2) then
            allocate(DevFuns%nx_all(size(funs%ufenlx_all,1),size(funs%ufenlx_all,2),size(funs%ufenlx_all,3)))
        else
            allocate(DevFuns%nx_all(size(funs%cvfenlx_all,1),size(funs%cvfenlx_all,2),size(funs%cvfenlx_all,3)))
        end if

    end subroutine allocate_multi_dev_shape_funs1

    subroutine allocate_multi_dev_shape_funs2(Mdims, GIdims, DevFuns, nx_all_FE_size)
        implicit none
        type (multi_dimensions), intent(in)  ::Mdims
        type(multi_GI_dimensions), intent(in) :: GIdims
        type (multi_dev_shape_funs), intent(inout) :: DevFuns
        logical, optional, intent(in) :: nx_all_FE_size!If true then the size of the generic nx_all is set like
                                                       !like ufenlx_all, otherwise, like cvfenlx_all
        !Local variable
        logical :: nx_all_FE_size2
        nx_all_FE_size2 = .false.

        allocate(DevFuns%detwei(GIdims%cv_ngi))
        allocate(DevFuns%ra(GIdims%cv_ngi))
        allocate(DevFuns%cvfenx_all(Mdims%Ndim, Mdims%cv_nloc,GIdims%cv_ngi))
        allocate(DevFuns%ufenx_all(Mdims%Ndim, Mdims%u_nloc, GIdims%cv_ngi))
        allocate(DevFuns%inv_jac(Mdims%Ndim, Mdims%Ndim, GIdims%cv_ngi))

        if (present(nx_all_FE_size)) then
            if (nx_all_FE_size) nx_all_FE_size2 = .true.
        end if

        if (nx_all_FE_size2) then
            allocate(DevFuns%nx_all(Mdims%Ndim, Mdims%u_nloc, GIdims%cv_ngi))
        else
            allocate(DevFuns%nx_all(Mdims%Ndim, Mdims%cv_nloc,GIdims%cv_ngi))
        end if

    end subroutine allocate_multi_dev_shape_funs2

    subroutine allocate_multi_dev_shape_funs3(cvfenlx_all, ufenlx_all, DevFuns, nx_all_FE_size)
        implicit none
        real, dimension(:,:,:), intent(in) :: cvfenlx_all, ufenlx_all
        type (multi_dev_shape_funs), intent(inout) :: DevFuns
        logical, optional, intent(in) :: nx_all_FE_size!If true then the size of the generic nx_all is set like
                                                       !like ufenlx_all, otherwise, like cvfenlx_all
        !Local variable
        logical :: nx_all_FE_size2
        nx_all_FE_size2 = .false.

        allocate(DevFuns%detwei(size(cvfenlx_all,3)))
        allocate(DevFuns%ra(size(cvfenlx_all,3)))
        allocate(DevFuns%cvfenx_all(size(cvfenlx_all,1),size(cvfenlx_all,2),size(cvfenlx_all,3)))
        allocate(DevFuns%ufenx_all(size(ufenlx_all,1),size(ufenlx_all,2),size(ufenlx_all,3)))
        allocate(DevFuns%inv_jac(size(ufenlx_all,1), size(ufenlx_all,1), size(ufenlx_all,3)))

        if (present(nx_all_FE_size)) then
            if (nx_all_FE_size) nx_all_FE_size2 = .true.
        end if

        if (nx_all_FE_size2) then
            allocate(DevFuns%nx_all(size(ufenlx_all,1),size(ufenlx_all,2),size(ufenlx_all,3))) ; DevFuns%nx_all=0.0
        else
            allocate(DevFuns%nx_all(size(cvfenlx_all,1),size(cvfenlx_all,2),size(cvfenlx_all,3))) ; DevFuns%nx_all=0.0
        end if
    end subroutine allocate_multi_dev_shape_funs3

    subroutine deallocate_multi_dev_shape_funs(DevFuns)
        implicit none
        type (multi_dev_shape_funs), intent(inout) :: DevFuns
        !Deallocate memory
        if (associated(DevFuns%detwei)) deallocate(DevFuns%detwei)
        if (associated(DevFuns%ra)) deallocate(DevFuns%ra)
        if (associated(DevFuns%ufenx_all)) deallocate(DevFuns%ufenx_all)
        if (associated(DevFuns%cvfenx_all)) deallocate(DevFuns%cvfenx_all)
        if (associated(DevFuns%nx_all)) deallocate(DevFuns%nx_all)
        if (associated(DevFuns%inv_jac)) deallocate(DevFuns%inv_jac)

        !Nullify pointers
        nullify(DevFuns%cvfenx_all);nullify(DevFuns%ufenx_all)
        nullify(DevFuns%detwei); nullify(DevFuns%ra)
        nullify(DevFuns%nx_all); nullify(DevFuns%inv_jac)
    end subroutine deallocate_multi_dev_shape_funs

    subroutine allocate_porous_adv_coefs(Mdims, upwnd)
        type (porous_adv_coefs), intent(inout) :: upwnd
        type (multi_dimensions), intent(in)  ::Mdims

!        if (.not.associated(upwnd%adv_coef)) allocate(upwnd%adv_coef(Mdims%ndim,Mdims%ndim,Mdims%nphase,Mdims%mat_nonods))
        if (.not.associated(upwnd%inv_adv_coef)) allocate(upwnd%inv_adv_coef(Mdims%ndim,Mdims%ndim,Mdims%nphase,Mdims%mat_nonods))
        if (.not.associated(upwnd%adv_coef_grad)) allocate(upwnd%adv_coef_grad(Mdims%ndim,Mdims%ndim,Mdims%nphase,Mdims%mat_nonods))
    end subroutine allocate_porous_adv_coefs

    subroutine deallocate_porous_adv_coefs(upwnd)
        type (porous_adv_coefs), intent(inout) :: upwnd

!        if (associated(upwnd%adv_coef)) deallocate(upwnd%adv_coef)!This memory is pointing to
                                            !multi_absorption%porousMedia as is being deallocated there
        if (associated(upwnd%inv_adv_coef)) deallocate(upwnd%inv_adv_coef)
        if (associated(upwnd%adv_coef_grad)) deallocate(upwnd%adv_coef_grad)

        nullify(upwnd%adv_coef); nullify(upwnd%inv_adv_coef);nullify(upwnd%adv_coef_grad)
    end subroutine deallocate_porous_adv_coefs


    subroutine allocate_multi_pipe_package(pipes, Mdims, Mspars)
        type (multi_pipe_package), intent(inout) :: pipes
        type (multi_dimensions), intent(in)  ::Mdims
        type (multi_sparsities), intent(in) :: Mspars

        if (Mdims%npres > 1) then
            if (.not.associated(pipes%gamma_pres_abs))        allocate( pipes%gamma_pres_abs( mdims%nphase,mdims%nphase,mdims%cv_nonods ))
            if (.not.associated(pipes%gamma_pres_abs_nano))   allocate( pipes%gamma_pres_abs_nano( mdims%nphase,mdims%nphase,mdims%cv_nonods ))
            if (.not.associated(pipes%mass_pipe))             allocate( pipes%mass_pipe( mdims%cv_nonods ))
            if (.not.associated(pipes%mass_cvfem2pipe))       allocate(pipes%mass_cvfem2pipe( mspars%cmc%ncol ))
            if (.not.associated(pipes%mass_pipe2cvfem))       allocate( pipes%mass_pipe2cvfem( mspars%cmc%ncol ))
            if (.not.associated(pipes%mass_cvfem2pipe_true))  allocate( pipes%mass_cvfem2pipe_true( mspars%cmc%ncol ))
        else
!            if (.not.associated(pipes%gamma_pres_abs))        allocate( pipes%gamma_pres_abs( 0,0,0 ))
!            if (.not.associated(pipes%gamma_pres_abs_nano))   allocate( pipes%gamma_pres_abs_nano( 0,0,0 ))
!            if (.not.associated(pipes%mass_pipe))             allocate( pipes%mass_pipe( 0 ))
!            if (.not.associated(pipes%mass_cvfem2pipe))       allocate(pipes%mass_cvfem2pipe( 0 ))
!            if (.not.associated(pipes%mass_pipe2cvfem))       allocate( pipes%mass_pipe2cvfem( 0 ))
!            if (.not.associated(pipes%mass_cvfem2pipe_true))  allocate( pipes%mass_cvfem2pipe_true( 0 ))
        end if
    end subroutine allocate_multi_pipe_package

    subroutine deallocate_multi_pipe_package(pipes)
        type (multi_pipe_package), intent(inout) :: pipes

        if (associated(pipes%gamma_pres_abs)) then
            deallocate( pipes%gamma_pres_abs); nullify(pipes%gamma_pres_abs)
        end if
        if (associated(pipes%gamma_pres_abs_nano)) then
            deallocate( pipes%gamma_pres_abs_nano); nullify(pipes%gamma_pres_abs_nano)
        end if
        if (associated(pipes%mass_pipe)) then
            deallocate(pipes%mass_pipe); nullify(pipes%mass_pipe)
        end if
        if (associated(pipes%mass_cvfem2pipe)) then
            deallocate(pipes%mass_cvfem2pipe); nullify(pipes%mass_cvfem2pipe)
        end if
        if (associated(pipes%mass_pipe2cvfem)) then
            deallocate(pipes%mass_pipe2cvfem); nullify(pipes%mass_pipe2cvfem)
        end if
        if (associated(pipes%mass_cvfem2pipe_true)) then
            deallocate(pipes%mass_cvfem2pipe_true); nullify(pipes%mass_cvfem2pipe_true)
        end if

    end subroutine deallocate_multi_pipe_package

    subroutine initialize_multi_outfluxes(outfluxes)
        implicit none
        type (multi_outfluxes), intent(inout) :: outfluxes
        !Local variables
        integer, dimension(2) :: shapes

        outfluxes%calculate_flux = have_option( "/io/dump_boundaryflux/surface_ids")
        ! Read in the surface IDs of the boundaries (if any) that you wish to integrate over into the (integer vector) variable outfluxes%outlet_id.
        ! No need to explicitly allocate outfluxes%outlet_id (done here internally)
        if (outfluxes%calculate_flux .and..not.(allocated(outfluxes%outlet_id))) then
            shapes = option_shape("/io/dump_boundaryflux/surface_ids")
            assert(shapes(1) >= 0)
            allocate(outfluxes%outlet_id(shapes(1)))
            call get_option( "/io/dump_boundaryflux/surface_ids", outfluxes%outlet_id)
        endif
        !At least size 1 to be used to calculate the whole mass of the domain, and keep valgrind happy!
        if (.not. allocated(outfluxes%outlet_id)) then
            allocate(outfluxes%outlet_id(1))
        end if

    end subroutine initialize_multi_outfluxes

    subroutine allocate_multi_outfluxes(Mdims, outfluxes)
        implicit none
        type (multi_dimensions), intent(in)  ::Mdims
        type (multi_outfluxes), intent(inout) :: outfluxes
        !Local variables
        integer :: k

        allocate(outfluxes%intflux(Mdims%nphase,size(outfluxes%outlet_id)))
        k = 1
        if (has_temperature) k = k + 1
        !(field -saturation, temperature-, Mdims%nphase, size(outfluxes%outlet_id))
        allocate(outfluxes%totout(k, Mdims%nphase, size(outfluxes%outlet_id)))

        outfluxes%intflux= 0.
        outfluxes%totout= 0.

    end subroutine allocate_multi_outfluxes

    subroutine destroy_multi_outfluxes(outfluxes)
        implicit none
        type (multi_outfluxes), intent(inout) :: outfluxes

        deallocate(outfluxes%totout, outfluxes%intflux,outfluxes%outlet_id )

    end subroutine destroy_multi_outfluxes

end module multi_data_types


