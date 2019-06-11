
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

module multi_magma

    use fldebug
    use futils
    use spud
    use fields
    use global_parameters, only: OPTION_PATH_LEN, PYTHON_FUNC_LEN, PI, is_porous_media
    use vector_tools

contains
    !========================================================
    !Functions Converting between Dimensional and Non-Dimensional Temperatures
    !========================================================
    real function temp_dim_to_nondim(T,Ts,Tl)
    implicit none
    real, intent(in):: T,Ts,Tl

        temp_dim_to_nondim = (T - Ts)/(Tl-Ts)

    end function
    !========================================================
    real function temp_nondim_to_dim(TPrime,Ts,Tl)
    implicit none
    real, intent(in):: TPrime,Ts,Tl

        temp_nondim_to_dim = TPrime*(Tl-Ts) + Ts

    end function
    !========================================================
    real function enth_dim_to_nondim(H,Hs,Hl)
    implicit none
    real, intent(in):: H,Hs,Hl

        enth_dim_to_nondim = (H - Hs)/(Hl-Hs)

    end function
    !========================================================
    real function enth_nondim_to_dim(HPrime,Hs,Hl)
    implicit none
    real, intent(in) :: HPrime,Hs,Hl

        enth_nondim_to_dim = HPrime*(Hl-Hs) + Hs

    end function

end module multi_magma
