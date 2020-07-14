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

SHELL = /bin/sh
PROGRAM = ../bin/b3d

FC      = mpif90
FCFLAGS = -I../include  -I/usr/lib/petscdir/petsc-3.6.4/include -I/usr/lib/petscdir/petsc-3.6.4/arch-linux2-c-debug/include -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/opal/mca/event/libevent2022/libevent -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/opal/mca/event/libevent2022/libevent/include -I/usr/lib/x86_64-linux-gnu/openmpi/include -ffast-math -frecord-marker=4  -I/usr/include -DHAVE_PYTHON -ffree-line-length-none -ffixed-line-length-none  -O3 -I/usr/lib/petscdir/petsc-3.6.4/include -I/usr/lib/petscdir/petsc-3.6.4/arch-linux2-c-debug/include -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/opal/mca/event/libevent2022/libevent -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/opal/mca/event/libevent2022/libevent/include -I/usr/lib/x86_64-linux-gnu/openmpi/include -fdefault-real-8 -I/usr/local/include -I./

MPIF90  = mpif90

LINKER  = mpicxx   

CC      = mpicc
CFLAGS  =  -I/usr/lib/petscdir/petsc-3.6.4/include -I/usr/lib/petscdir/petsc-3.6.4/arch-linux2-c-debug/include -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/opal/mca/event/libevent2022/libevent -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/opal/mca/event/libevent2022/libevent/include -I/usr/lib/x86_64-linux-gnu/openmpi/include -g -O2 -DHAVE_PYTHON -I/usr/include/python2.7 -I../include -I/usr/local/include -g

CXX     = mpicxx
CXXFLAGS=  -I/usr/lib/petscdir/petsc-3.6.4/include -I/usr/lib/petscdir/petsc-3.6.4/arch-linux2-c-debug/include -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/opal/mca/event/libevent2022/libevent -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/opal/mca/event/libevent2022/libevent/include -I/usr/lib/x86_64-linux-gnu/openmpi/include -g -O2 -DHAVE_PYTHON -I/usr/include/python2.7 -I../include -I/usr/local/include

MAKE    = make
AR  = ar
ARFLAGS = cr

LIB = ../lib/libb3d.a

OBJS    =  frame.o  GID_B3D.o 

.SUFFIXES: .f90 .F90 .c .cpp .o .a

.f90.o:
	@echo "        FC $<"
	$(FC) $(FCFLAGS) -c $< || (rm $(LIB) && false)
.F90.o:
	@echo "        FC $<"
	$(FC) $(FCFLAGS) -c $< || (rm $(LIB) && false)
.c.o:
	@echo "        CC $<"
	$(CC) $(CFLAGS) -c $< || (rm $(LIB) && false)
.cpp.o:
	@echo "        CXX $<"
	$(CXX) $(CXXFLAGS) -c $< || (rm $(LIB) && false)

default: $(OBJS)
	$(AR) $(ARFLAGS) $(LIB) $(OBJS)
	@$(LINKER)  -I/usr/lib/petscdir/petsc-3.6.4/include -I/usr/lib/petscdir/petsc-3.6.4/arch-linux2-c-debug/include -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/opal/mca/event/libevent2022/libevent -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/opal/mca/event/libevent2022/libevent/include -I/usr/lib/x86_64-linux-gnu/openmpi/include -I./include -o $(PROGRAM) $(LIB)

# cancel implicit rule that tries to make .o from .mod with m2c:
%.o: %.mod

clean:
	rm -f *.o *.d 

