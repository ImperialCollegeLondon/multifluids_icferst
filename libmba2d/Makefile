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
#    C.Pain@Imperial.ac.uk
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

FC      = mpif90
FCFLAGS = -I../include -ffast-math -frecord-marker=4 -I/usr/include -DHAVE_NUMPY -I/usr/lib/python3/dist-packages/numpy/core/include -fno-realloc-lhs -ffree-line-length-none -ffixed-line-length-none  -O3 -I/usr/lib/petscdir/3.19/include -I/usr/lib/petscdir/3.19//include -I/usr/include/python3.12 -fdefault-real-8 -fdefault-double-8

F77      = mpif90
FFLAGS = -I../include -ffast-math -frecord-marker=4  -DHAVE_NUMPY -I/usr/lib/python3/dist-packages/numpy/core/include -fno-realloc-lhs -ffree-line-length-none -ffixed-line-length-none  -O3 -fdefault-real-8 -fdefault-double-8

UTFC    = ../bin/unittest_f90

CC      = mpicc
CFLAGS  =  -I/usr/include/python3.12 -I/usr/include/python3.12 -DHAVE_NUMPY -I/usr/lib/python3/dist-packages/numpy/core/include -DFLUIDITY_PYTHON_MODULE_PATH=\"/home/mbahlali/Softwares/MultiFluids_Dev_ubuntu24/python\"  -O3 -I/usr/lib/petscdir/3.19/include -I/usr/lib/petscdir/3.19//include

CXX     = mpicxx
CXXFLAGS= -std=c++11

MAKE    = gmake
AR  = ar
ARFLAGS = cr

LIB = ../lib/libmba2d.a

TESTLIBS = $(shell echo -lmba2d -lzoltan -lparmetis -lmetis -L/usr/lib/petscdir/3.19//lib -lpetsc -lnetcdff -lnetcdf  -llapack -ludunits2 -lpthread -lm -lstdc++  -L/usr/lib/python3.12/config-3.12-x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -lpython3.12 -ldl  -lm  -ldl -lpthread -lutil -lm  -L/usr/lib/gcc/x86_64-linux-gnu/13 -L/usr/lib/gcc/x86_64-linux-gnu/13/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/13/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/13/../../.. -lgfortran -lm -lquadmath  -L/usr/lib/gcc/x86_64-linux-gnu/13 -L/usr/lib/gcc/x86_64-linux-gnu/13/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/13/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/13/../../.. -lgfortran -lm -lquadmath -L./lib  -lvtkCommonCore-9.1 -lvtkCommonDataModel-9.1 -lvtkCommonExecutionModel-9.1 -lvtkFiltersCore-9.1 -lvtkFiltersVerdict-9.1 -lvtkIOCore-9.1 -lvtkIOXML-9.1 -lvtkIOGeometry-9.1 -lvtkIOLegacy-9.1 -lvtkIOParallelXML-9.1 -lvtkParallelMPI-9.1 -lvtkIOXMLParser-9.1 -lvtkCommonSystem-9.1 -lvtkCommonMisc-9.1 -lvtkCommonMath-9.1 -lvtksys-9.1 -ldl -lpthread -L/usr/lib/x86_64-linux-gnu/openmpi/lib -lmpi_cxx -lmpi  | sed 's@./lib/lib\([a-z]*\)\.a@-l\1@g')

# These are all the objects that are supposed to go into the makefile
OBJS = ani2.o auxSE.o auxSF.o auxSP.o calcrv.o clpsF.o delaunay.o dsort.o error.o forlibmba.o \
       insrtP.o lintrp2D.o list_bsct.o loadM.o makM.o makQ_linear.o makS.o mba_analytic.o \
       mba_nodal.o minim.o moveP.o nlnfnc.o prjcrv.o saveM.o smoothing.o splitE.o \
       swapF.o time.o untangle.o update.o utils.o ZZ.o isnan.o isnan_f.o

.SUFFIXES: .f90 .F90 .c .o .a .f

.f.o:
	@echo "        FC $<"
	$(FC) $(FCFLAGS) -c $<
.f90.o:
	@echo "        FC $<"
	$(FC) $(FCFLAGS) -c $<
.F90.o:
	@echo "        FC $<"
	$(FC) $(FCFLAGS) -c $<
.cpp.o:
	@echo "        CXX $<"
	$(CXX) $(CXXFLAGS) -I../include -c $<

default: $(OBJS)
	$(AR) $(ARFLAGS) $(LIB) $(OBJS)
	cp *.mod ../include/

clean:
	rm -f *.o *.mod

