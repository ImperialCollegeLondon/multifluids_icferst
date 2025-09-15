#    Copyright (C) 2006 Imperial College London and others.
#
#    Please see the AUTHORS file in the main source directory for a full list
#    of copyright holders.
#
#    Prof. C Pain
#    Applied Modelling and Computation Group
#    Department of Earth Science and Engineering
#    Imperial College London
#
#    amcgsoftware@imperial.ac.uk
#
#    This library is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation,
#    version 2.1 of the License.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with this library; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
#    USA
SHELL = /bin/bash

PACKAGE_NAME = fluidity
FLUIDITY = $(PACKAGE_NAME)

LIB = ../lib/lib$(FLUIDITY).a

AR = ar
ARFLAGS = cr

CC  = mpicc
CXX = mpicxx
FC = mpif90

LINKER = $(CXX)

# Adding /usr/include like this is not ideal, see:
# http://gcc.gnu.org/bugzilla/show_bug.cgi?id=35707
# If /usr/include occurs too early it's going to
# interfere with other compilers - and it would completely
# stuff C++ with obsecure, ungoogleable, error messages. 
FCFLAGS   = -I../include @MODFLAG@../include -ffast-math -frecord-marker=4  -DHAVE_NUMPY -I/usr/lib/python3/dist-packages/numpy/core/include -fno-realloc-lhs -ffree-line-length-none -ffixed-line-length-none  -O3 -I/media/hh210/data2/Petsc/petsc-3.19.5/include -I/media/hh210/data2/Petsc/petsc-3.19.5/arch-linux-c-debug/include -I/usr/include/python3.12 -fdefault-real-8 -fdefault-double-8  -I/usr/include/python3.12 -I/usr/include/python3.12 -DHAVE_NUMPY -I/usr/lib/python3/dist-packages/numpy/core/include -I/media/hh210/data2/Petsc/petsc-3.19.5/include -I/media/hh210/data2/Petsc/petsc-3.19.5/arch-linux-c-debug/include -DHAVE_PETSC  -I/usr/include/vtk-9.1 -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -I/usr/include

CXXFLAGS = -I../include -std=c++11  -I/usr/include/python3.12 -I/usr/include/python3.12 -DHAVE_NUMPY -I/usr/lib/python3/dist-packages/numpy/core/include -I/media/hh210/data2/Petsc/petsc-3.19.5/include -I/media/hh210/data2/Petsc/petsc-3.19.5/arch-linux-c-debug/include -DHAVE_PETSC  -I/usr/include/vtk-9.1 -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi

LDFLAGS   =   
LIBS    = -L../lib -l$(FLUIDITY) -lmba2d -lparmetis -lmetis -L/media/hh210/data2/Petsc/petsc-3.19.5/arch-linux-c-debug/lib -lpetsc -lnetcdff -lnetcdf  -llapack -ludunits2 -lpthread -lm -lstdc++  -L/usr/lib/python3.12/config-3.12-x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -lpython3.12 -ldl  -lm  -ldl -lpthread -lutil -lm  -L/usr/lib/gcc/x86_64-linux-gnu/13 -L/usr/lib/gcc/x86_64-linux-gnu/13/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/13/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/13/../../.. -lgfortran -lm -lquadmath  -L/usr/lib/gcc/x86_64-linux-gnu/13 -L/usr/lib/gcc/x86_64-linux-gnu/13/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/13/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/13/../../.. -lgfortran -lm -lquadmath -L./lib  -lvtkCommonCore-9.1 -lvtkCommonDataModel-9.1 -lvtkCommonExecutionModel-9.1 -lvtkFiltersCore-9.1 -lvtkFiltersVerdict-9.1 -lvtkIOCore-9.1 -lvtkIOXML-9.1 -lvtkIOGeometry-9.1 -lvtkIOLegacy-9.1 -lvtkIOParallelXML-9.1 -lvtkParallelMPI-9.1 -lvtkIOXMLParser-9.1 -lvtkCommonSystem-9.1 -lvtkCommonMisc-9.1 -lvtkCommonMath-9.1 -lvtksys-9.1 -ldl -lpthread -L/usr/lib/x86_64-linux-gnu/openmpi/lib -lmpi_cxx -lmpi  -lblas  -L/usr/lib/gcc/x86_64-linux-gnu/13 -L/usr/lib/gcc/x86_64-linux-gnu/13/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/13/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/13/../../.. -lgfortran -lm -lquadmath \
		  ../lib/libvtkfortran.a /media/hh210/data2/ICFEMRES/ICFERST_petsc3.19/lib/libspatialindex.a /media/hh210/data2/ICFEMRES/ICFERST_petsc3.19/lib/libspud.a ./lib/libjudy.a

OBJS = NetCDF_reader.o SampleNetCDF2.o import_bath_data.o read_netcdf_interface.o

.SUFFIXES: .c .o .cpp .F90 .a .so

.F90.o:
	@echo "        FC $<"
	$(FC) $(FCFLAGS) $(GENFLAGS) -c $< 
.c.o:
	@echo "        CC $<"
	$(CC) $(CFLAGS) -c $<
.cpp.o:
	@echo "        CXX $<"
	$(CXX) $(CXXFLAGS) -c $<

$(LIB): $(OBJS)
	@mkdir -p ../lib

default: $(LIB)

DEPTH_OBJS=import_bath_data.cpp

install:

clean:
	rm -f *.o *.mod
include Makefile.dependencies
