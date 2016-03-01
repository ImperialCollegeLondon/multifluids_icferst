  
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

    use global_parameters, only: option_path_len, is_porous_media

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
        integer :: mat_nloc   !??
        integer :: totele     !Total number of elements
        integer :: stotel     !Total number of surface elements?
        integer :: cv_nonods  !Total number of control volumes
        integer :: p_nonods   !Total number of pressure nodes
        integer :: mat_nonods !Total number of ???
        integer :: u_nonods   !Total number of velocity nodes
        integer :: xu_nonods  !Total number of velocity nodes of the Continuous mesh
        integer :: x_nonods   !Total number of control volumes of the Continuous mesh
        integer :: ph_nloc    !Number of ????
        integer :: ph_nonods  !Total number of ????
        integer :: nphase     !Total number of phases
        integer :: npres      !Total number of pressure
        integer :: n_in_pres  !nphase/npres

    end type multi_dimensions

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
    end type multi_shape_funs

    !This type comprises the four necessary variables to represent matrices using a CSR structure
    type multi_sparsity
        integer :: ncol
        integer, pointer, dimension(:) :: fin=> null()
        integer, pointer, dimension(:) :: col=> null()
        integer, pointer, dimension(:) :: mid=> null()
    end type multi_sparsity

    !This data type contains all the sparcities necessary in the multiphase prototype code
    type multi_sparsities
        type (multi_sparsity) :: acv     !CV multi-phase eqns (e.g. vol frac, temp)
        type (multi_sparsity) :: small_acv !Local CV multi-phase eqns (e.g. vol frac, temp)
        type (multi_sparsity) :: mcy     !Force balance plus cty multi-phase eqns
        type (multi_sparsity) :: ele     !Element connectivity
        type (multi_sparsity) :: dgm_pha !Force balance sparsity
        type (multi_sparsity) :: ct      !CT sparsity - global continuity eqn
        type (multi_sparsity) :: C       !C sparsity operating on pressure in force balance
        type (multi_sparsity) :: cmc     !pressure matrix for projection method
        type (multi_sparsity) :: m       !CV-FEM matrix
        type (multi_sparsity) :: ph      !ph matrix
    end type multi_sparsities


contains
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
    end subroutine deallocate_multi_shape_funs


    !This subroutine, despite it can be called by itself it is highly recommended to be called ONLY through multi_sparsity/Get_Sparsity_Patterns
    subroutine allocate_multi_sparsities(Mspars, Mdims, mx_ncolacv, mx_ncolmcy, nlenmcy, mx_ncoldgm_pha, mx_nct, mx_nc, mx_ncolm, mx_ncolph)
        !This subroutine allocates part of the memory inside Mspars
        implicit none
        type (multi_sparsities), intent(inout) :: Mspars
        type(multi_dimensions), intent(in) :: Mdims
        integer :: mx_ncolacv, mx_ncolmcy, nlenmcy, mx_ncoldgm_pha, mx_nct, mx_nc, mx_ncolm, mx_ncolph

        allocate( Mspars%ACV%fin( Mdims%cv_nonods * Mdims%nphase + 1 ), Mspars%ACV%col( mx_ncolacv ), Mspars%ACV%mid( Mdims%cv_nonods * Mdims%nphase ), &
            Mspars%MCY%fin( nlenmcy + 1 ), Mspars%MCY%col( mx_ncolmcy ), Mspars%MCY%mid( nlenmcy ), &
            Mspars%DGM_PHA%fin( Mdims%u_nonods * Mdims%nphase * Mdims%ndim + 1 ), Mspars%DGM_PHA%col( mx_ncoldgm_pha ), &
            Mspars%DGM_PHA%mid( Mdims%u_nonods * Mdims%nphase * Mdims%ndim ), &
            Mspars%CT%fin( Mdims%cv_nonods + 1 ), Mspars%CT%col( mx_nct ), &
            Mspars%C%fin( Mdims%u_nonods + 1 ), Mspars%C%col( mx_nc ), &
            Mspars%CMC%fin( Mdims%cv_nonods + 1 ), Mspars%CMC%col( 0 ), Mspars%CMC%mid( Mdims%cv_nonods ), &
            Mspars%M%fin( Mdims%cv_nonods + 1 ), Mspars%M%col( mx_ncolm ), Mspars%M%mid( Mdims%cv_nonods ), &
            Mspars%ph%fin( Mdims%ph_nonods + 1 ), Mspars%ph%col( mx_ncolph ) )

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

        if (associated(Mspars%mcy%fin))       deallocate(Mspars%mcy%fin)
        if (associated(Mspars%mcy%col))       deallocate(Mspars%mcy%col)
        if (associated(Mspars%mcy%mid))       deallocate(Mspars%mcy%mid)

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

!        deallocate(Mspars%ele%fin, Mspars%ele%col, Mspars%ele%mid)!Not deallocate this one since it points to state
        !Proceed to nullify
        nullify(Mspars%acv%fin, Mspars%acv%col, Mspars%acv%mid)
        nullify(Mspars%small_acv%fin, Mspars%small_acv%col, Mspars%small_acv%mid)
        nullify(Mspars%mcy%fin, Mspars%mcy%col, Mspars%mcy%mid)
        nullify(Mspars%ele%fin, Mspars%ele%col, Mspars%ele%mid)
        nullify(Mspars%dgm_pha%fin, Mspars%dgm_pha%col, Mspars%dgm_pha%mid)
        nullify(Mspars%ct%fin, Mspars%ct%col, Mspars%ct%mid)
        nullify(Mspars%C%fin, Mspars%C%col, Mspars%C%mid)
        nullify(Mspars%CMC%fin, Mspars%CMC%col, Mspars%CMC%mid)
        nullify(Mspars%M%fin, Mspars%M%col, Mspars%M%mid)
        nullify(Mspars%ph%fin, Mspars%ph%col, Mspars%ph%mid)


    end subroutine deallocate_multi_sparsities


end module multi_data_types


