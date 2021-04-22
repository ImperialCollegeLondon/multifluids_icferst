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

FC      = mpif90
FCFLAGS = -I../include -J../include -ffast-math -frecord-marker=4  -fno-realloc-lhs -ffree-line-length-none -ffixed-line-length-none  -O3 -I/usr/lib/petscdir/3.8.3/include -I/usr/lib/petscdir/3.8.3/linux-gnu-c-opt/include -I/usr/include/hdf5/openmpi -I/home/amin/Desktop/IC-FIRST/include/ -fdefault-real-8 -fdefault-double-8 -I/usr/local/include -I./

CC      = mpicc
CFLAGS  =  -I/usr/include/python2.7  -O3 -I../include -I/usr/local/include

CXX     = mpicxx
CXXFLAGS=  -I/usr/include/python2.7  -O3 -I../include -I/usr/local/include

MAKE    = make
AR  = ar
ARFLAGS = cr

LIB = ../lib/lib$(FLUIDITY).a

OBJS = DQMOM.o \

.SUFFIXES: .f90 .F90 .c .o .a

.f90.o:
	@echo "        FC $<"
	$(FC) $(FCFLAGS) -c $<
.F90.o:
	@echo "        FC $<"
	$(FC) $(FCFLAGS) $(GENFLAGS) -c $< 

$(LIB): $(OBJS)
	@mkdir -p ../lib

default: $(LIB)

clean:
	rm -f *.o *.d *.mod

include Makefile.dependencies
