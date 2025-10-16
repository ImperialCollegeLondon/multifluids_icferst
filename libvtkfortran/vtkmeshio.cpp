/* Copyright (C) 2004- Imperial College London and others.
   
   Please see the AUTHORS file in the main source directory for a full
   list of copyright holders.
   
   Adrian Umpleby
   Applied Modelling and Computation Group
   Department of Earth Science and Engineering
   Imperial College London
   
   adrian@imperial.ac.uk
   
   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License.
   
   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.
   
   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
   USA
*/

/* Copyright (C) 2004- Imperial College London and others.
   
   Please see the AUTHORS file in the main source directory for a full
   list of copyright holders.
   
   Adrian Umpleby
   Applied Modelling and Computation Group
   Department of Earth Science and Engineering
   Imperial College London
   
   adrian@imperial.ac.uk
   
   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License.
   
   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.
   
   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
   USA
*/

#include "confdefs.h"

#ifdef HAVE_VTK
#include <vtkVersion.h>
#if VTK_VERSION_NUMBER >= 89000000000
#include "meshio/vtk9meshio.cpp"
#else
#include "meshio/vtk7meshio.cpp"
#endif
#else
#include "meshio/vtkmeshio-dummy.cpp"
#endif
