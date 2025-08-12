# Copyright (C) 2004-2006 by Gerard Gorman
# Copyright (C) 2006- Imperial College London and others.
#
#  Please see the AUTHORS file in the main source directory for a full list
#  of copyright holders.
#   
#  Dr Gerard J Gorman
#  Applied Modelling and Computation Group
#  Department of Earth Science and Engineering
#  Imperial College London
#   
#  g.gorman@imperial.ac.uk
#   
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License.
#   
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#  
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
#  USA

SHELL = /bin/bash

FC = mpif90
FCFLAGS = -I./include  -I/usr/include/python3.12 -I/usr/include/python3.12 -DHAVE_NUMPY -I/usr/lib/python3/dist-packages/numpy/core/include -I/usr/lib/petscdir/3.19/include -I/usr/lib/petscdir/3.19//include -DHAVE_PETSC  -I/usr/include/vtk-9.1 -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -ffast-math -frecord-marker=4 -I/usr/include -DHAVE_NUMPY -I/usr/lib/python3/dist-packages/numpy/core/include -fno-realloc-lhs -ffree-line-length-none -ffixed-line-length-none  -O3 -I/usr/lib/petscdir/3.19/include -I/usr/lib/petscdir/3.19//include -I/usr/include/python3.12 -fdefault-real-8 -fdefault-double-8

CXX = mpicxx
CXXFLAGS = -I./include  -I/usr/include/python3.12 -I/usr/include/python3.12 -DHAVE_NUMPY -I/usr/lib/python3/dist-packages/numpy/core/include -I/usr/lib/petscdir/3.19/include -I/usr/lib/petscdir/3.19//include -DHAVE_PETSC  -I/usr/include/vtk-9.1 -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -std=c++11

LINKER = mpicxx
LIBS = -lmba2d -lzoltan -lparmetis -lmetis -L/usr/lib/petscdir/3.19//lib -lpetsc -lnetcdff -lnetcdf  -llapack -ludunits2 -lpthread -lm -lstdc++  -L/usr/lib/python3.12/config-3.12-x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -lpython3.12 -ldl  -lm  -ldl -lpthread -lutil -lm  -L/usr/lib/gcc/x86_64-linux-gnu/13 -L/usr/lib/gcc/x86_64-linux-gnu/13/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/13/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/13/../../.. -lgfortran -lm -lquadmath  -L/usr/lib/gcc/x86_64-linux-gnu/13 -L/usr/lib/gcc/x86_64-linux-gnu/13/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/13/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/13/../../.. -lgfortran -lm -lquadmath -L./lib  -lvtkCommonCore-9.1 -lvtkCommonDataModel-9.1 -lvtkCommonExecutionModel-9.1 -lvtkFiltersCore-9.1 -lvtkFiltersVerdict-9.1 -lvtkIOCore-9.1 -lvtkIOXML-9.1 -lvtkIOGeometry-9.1 -lvtkIOLegacy-9.1 -lvtkIOParallelXML-9.1 -lvtkParallelMPI-9.1 -lvtkIOXMLParser-9.1 -lvtkCommonSystem-9.1 -lvtkCommonMisc-9.1 -lvtkCommonMath-9.1 -lvtksys-9.1 -ldl -lpthread -L/usr/lib/x86_64-linux-gnu/openmpi/lib -lmpi_cxx -lmpi 

AR = ar
ARFLAGS = cr

MAKE = gmake

LIB = libvtkfortran.a
SLIB = libvtkfortran.so
DLIB = libvtkfortran.dylib

MAJOR = 1
VERSION = 1

OBJS = vtkfortran.o vtkmeshio.o fvtkfortran.o
TINY_OBJS = tinyxmlparser.o tinyxmlerror.o tinyxml.o

.SUFFIXES:
.SUFFIXES: .F90 .cpp .o .a

.cpp.o:
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $<

.F90.o:
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $<

default: $(OBJS)
	$(AR) $(ARFLAGS) $(LIB) $(OBJS)

shared: $(OBJS) $(TINY_OBJS)

	$(LINKER) -shared -Wl,-soname,$(SLIB).$(MAJOR),-z,defs -o $(SLIB).$(VERSION) $(OBJS) $(TINY_OBJS) $(LIBS)

dynamic: $(OBJS)
	$(CXX) -dynamiclib -install_name $(DLIB) -o $(DLIB) $(OBJS) $(LIBS)

clean:
	rm -f *.o *.mod $(LIB) $(SLIB)* $(DLIB)

distclean: clean
	rm -f make.log include/confdefs.h include/vtk.h
