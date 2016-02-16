  
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


    !Data structure to store all the shape functions to facilitate its movement throughtout the code
    type multi_shape_funs
        real, allocatable, dimension( : , : ) :: cvn ! dimension( cv_nloc, cv_ngi )
        real, allocatable, dimension( : ) :: cvweight!dimension( cv_ngi )
        real, allocatable, dimension(  : , : ) :: cvfen!dimension( cv_nloc, cv_ngi )
        real, allocatable, dimension(  : , : , :)  ::  cvfenlx_all!dimension( ndim, cv_nloc, cv_ngi )
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
    end type multi_shape_funs


contains

    subroutine allocate_multi_shape_funs(shape_fun, ndim, cv_ngi, cv_nloc, u_nloc, nface, scvngi, sbcvngi, cv_snloc, u_snloc, ncolgpts )
    !This subroutine allocates all the arrays in a multi_shape_funs data type
        implicit none
        integer, intent(in) :: ndim, cv_ngi, cv_nloc, u_nloc, nface, scvngi, sbcvngi, cv_snloc, u_snloc, ncolgpts
        type(multi_shape_funs), intent(inout) :: shape_fun

        !Proceed to allocate the variables
        allocate( shape_fun%cvn(cv_nloc, cv_ngi) )
        allocate(shape_fun%cvweight(cv_ngi))
        allocate(shape_fun%cvfen(cv_nloc, cv_ngi ))
        allocate(shape_fun%cvfenlx_all( ndim, cv_nloc, cv_ngi ))
        allocate(shape_fun%ufen( u_nloc, cv_ngi ))
        allocate(shape_fun%ufen( u_nloc, cv_ngi ))
        allocate(shape_fun%ufenlx_all( ndim, u_nloc, cv_ngi ))
        allocate(shape_fun%cv_neiloc( cv_nloc, scvngi ))
        allocate(shape_fun%cv_on_face( cv_nloc, scvngi ))
        allocate(shape_fun%cvfem_on_face( cv_nloc, scvngi ))
        allocate(shape_fun%scvfen( cv_nloc, scvngi ))
        allocate(shape_fun%scvfenslx( cv_nloc, scvngi ))
        allocate(shape_fun%scvfensly( cv_nloc, scvngi ))
        allocate(shape_fun%scvfeweigh( scvngi ))
        allocate(shape_fun%scvfenlx_all(ndim, cv_nloc, scvngi))
        allocate(shape_fun%sufen( u_nloc, scvngi ))
        allocate(shape_fun%sufenslx( u_nloc, scvngi ))
        allocate(shape_fun%sufensly( u_nloc, scvngi ))
        allocate(shape_fun%sufenlx_all( ndim, u_nloc, scvngi ))
        allocate(shape_fun%u_on_face( u_nloc, scvngi ))
        allocate(shape_fun%ufem_on_face( u_nloc, scvngi ))
        allocate(shape_fun%sbcvn(cv_snloc, sbcvngi))
        allocate(shape_fun%sbcvfen( cv_snloc, sbcvngi ))
        allocate(shape_fun%sbcvfenslx( cv_snloc, sbcvngi ))
        allocate(shape_fun%sbcvfensly( cv_snloc, sbcvngi ))
        allocate(shape_fun%sbcvfeweigh(sbcvngi))
        allocate(shape_fun%sbcvfenlx_all( ndim, cv_snloc, sbcvngi ))
        allocate(shape_fun%sbufen(u_snloc, sbcvngi))
        allocate(shape_fun%sbufenslx(u_snloc, sbcvngi))
        allocate(shape_fun%sbufensly(u_snloc, sbcvngi))
        allocate(shape_fun%sbufenlx_all( ndim, u_snloc, sbcvngi ))
        allocate(shape_fun%cv_sloclist( nface, cv_snloc ))
        allocate(shape_fun%u_sloclist( nface, u_snloc ))
        allocate(shape_fun%findgpts(cv_nloc + 1))
        allocate(shape_fun%colgpts( cv_nloc * scvngi ))


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


