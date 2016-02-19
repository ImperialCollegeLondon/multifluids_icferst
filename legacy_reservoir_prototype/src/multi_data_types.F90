  
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
        real, allocatable, dimension( : , : ) :: cvn ! dimension( cv_nloc, cv_ngi )
        real, allocatable, dimension( : ) :: cvweight!dimension( cv_ngi )
        real, allocatable, dimension(  : , : ) :: cvfen!dimension( cv_nloc, cv_ngi )
        real, allocatable, dimension( : , : , :)  ::  cvfenlx_all!dimension( ndim, cv_nloc, cv_ngi )
        real, allocatable, dimension(  : , : )  :: ufen!dimension( u_nloc, cv_ngi )
        real, allocatable, dimension(  : , :,: )  :: ufenlx_all!dimension( ndim, u_nloc, cv_ngi )
        integer, allocatable, dimension(  : , : )  :: cv_neiloc!dimension( cv_nloc, scvngi )
        logical, allocatable, dimension(  : , : ) :: cv_on_face, cvfem_on_face!dimension( cv_nloc, scvngi )
        real, allocatable, dimension(  : , : )  :: scvfen, scvfenslx, scvfensly!dimension( cv_nloc, scvngi )
        real, allocatable, dimension( : )  :: scvfeweigh!dimension( scvngi )
        real, allocatable, dimension(  : , : ,: )  :: scvfenlx_all!dimension( ndim, cv_nloc, scvngi )
        real, allocatable, dimension(  : , : )  :: sufen, sufenslx, sufensly!dimension( u_nloc, scvngi )
        real, allocatable, dimension(  : , :, : )  :: sufenlx_all !dimension( ndim, u_nloc, scvngi )
        logical, allocatable, dimension(  : , : )  :: u_on_face, ufem_on_face!dimension( u_nloc, scvngi )
        real, allocatable, dimension(  : , : )  :: sbcvn!dimension( cv_snloc, sbcvngi )
        real, allocatable, dimension(  : , : )  :: sbcvfen, sbcvfenslx, sbcvfensly! dimension( cv_snloc, sbcvngi )
        real, allocatable, dimension( : )  :: sbcvfeweigh!dimension( sbcvngi )
        real, allocatable, dimension(  : , :, : )  :: sbcvfenlx_all!dimension( ndim, cv_snloc, sbcvngi )
        real, allocatable, dimension(  : , : )  :: sbufen, sbufenslx, sbufensly!dimension( u_snloc, sbcvngi )
        real, allocatable, dimension(  : , :,: )  :: sbufenlx_all !dimension( ndim, u_snloc, sbcvngi )
        integer, allocatable, dimension(  : , : )  :: cv_sloclist!dimension( nface, cv_snloc )
        integer, allocatable, dimension(  : , : )  :: u_sloclist!dimension( nface, u_snloc )
        integer, allocatable, dimension( : )  :: findgpts!dimension( cv_nloc + 1 )
        integer, allocatable, dimension( : )  :: colgpts!dimension( cv_nloc * scvngi )
        integer :: ncolgpts
    end type multi_shape_funs


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
        deallocate( shape_fun%cvn)
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

    end subroutine deallocate_multi_shape_funs
end module multi_data_types


