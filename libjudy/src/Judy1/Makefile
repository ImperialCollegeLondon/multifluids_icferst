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

MPIF90  = mpif90

CC      = mpicc
CFLAGS  = -I/usr/include/vtk-5.8 -I/usr/include/python2.7 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include -I/usr/lib/petscdir/3.6.3/include -I/usr/lib/petscdir/3.6.3/linux-gnu-c-opt/include -I/usr/include -I/usr/lib/openmpi/include -I/usr/lib/openmpi/include/openmpi -Iinclude/ -DHAVE_PETSC  -I/usr/include/vtk-5.8 -DHAVE_VTK=1 -DHAVE_VTK=1 -fbounds-check -Wall -g -O0 -I/usr/include/python2.7 -DHAVE_NUMPY -I/usr/lib/python2.7/dist-packages/numpy/core/include  -DJU_64BIT -I../include -I/usr/local/include -I.. -I../JudyCommon -DJUDY1

MAKE    = make

LIBS    = -L../lib -l$(FLUIDITY) -lmba2d -lvtkIO -lvtkGraphics -lvtkFiltering -lvtkCommon -ldl -lzoltan -Wl,-rpath,/usr/lib/petscdir/3.6.3/linux-gnu-c-opt/lib -Wl,-rpath,/usr/lib/petscdir/3.6.3/linux-gnu-c-opt/lib -L/usr/lib/petscdir/3.6.3/linux-gnu-c-opt/lib -lpetsc -Wl,-rpath,/usr/lib/petscdir/3.6.3/linux-gnu-c-opt/lib -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lumfpack -lklu -lcholmod -lbtf -lccolamd -lcolamd -lcamd -lamd -lsuitesparseconfig -lml -lHYPRE -lscalapack -lexoIIv2for -lexodus -llapack -lblas -lparmetis -Wl,-rpath,/usr/lib -L/usr/lib -lnetcdf -ltriangle -lmetis -lchaco -lctetgen -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lX11 -lptesmumps -lptscotch -lptscotcherr -lscotch -lscotcherr -lm -Wl,-rpath,/usr/lib/openmpi/lib -L/usr/lib/openmpi/lib -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/4.8 -L/usr/lib/gcc/x86_64-linux-gnu/4.8 -Wl,-rpath,/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -lmpi_f90 -lmpi_f77 -lgfortran -lm -lgfortran -lm -lquadmath -lm -lmpi_cxx -lstdc++ -lrt -lm -lz -Wl,-rpath,/usr/lib/openmpi/lib -L/usr/lib/openmpi/lib -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/4.8 -L/usr/lib/gcc/x86_64-linux-gnu/4.8 -Wl,-rpath,/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -Wl,-rpath,/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -ldl -lmpi -lhwloc -lgcc_s -lpthread -ldl -lnetcdf -llapack -ludunits2 -lpthread -lm -lstdc++ -L/usr/lib -lpython2.7 -Xlinker -export-dynamic -Wl,-O0 -Wl,-Bsymbolic-functions -L/usr/lib -lz -lpthread -ldl -lutil -L/usr/lib/gcc/x86_64-linux-gnu/4.8 -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../.. -lgfortran -lm -lquadmath -L/usr/lib/gcc/x86_64-linux-gnu/4.8 -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../.. -lgfortran -lm -lquadmath -L./lib  -lblas  -L/usr/lib/gcc/x86_64-linux-gnu/4.8 -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../.. -lgfortran -lm -lquadmath

OBJS = Judy1Tables.o Judy1Test.o j__udy1Test.o Judy1Set.o Judy1SetArray.o Judy1Unset.o Judy1First.o Judy1Next.o \
			 Judy1Prev.o Judy1NextEmpty.o Judy1PrevEmpty.o Judy1Count.o Judy1ByCount.o Judy1FreeArray.o Judy1MemUsed.o \
			 Judy1MemActive.o Judy1Cascade.o Judy1Decascade.o Judy1CreateBranch.o Judy1InsertBranch.o Judy1MallocIF.o

.SUFFIXES: .f90 .F90 .cpp .c .o .a

# cancel implicit rule that tries to make .o from .mod with m2c:
%.o: %.mod

default: $(OBJS)
	cp *.o ../obj/

.c.o:
	@echo "        CC $<"
	$(CC) $(CFLAGS) -c $<

j__udy1Test.o: j__udy1Test.c
	@echo "        CC $<"
	$(CC) $(CFLAGS) -DJUDYGETINLINE -c $<

Judy1Next.o: Judy1Next.c
	@echo "        CC $<"
	$(CC) $(CFLAGS) -DJUDYNEXT -c $<

Judy1NextEmpty.o: Judy1NextEmpty.c
	@echo "        CC $<"
	$(CC) $(CFLAGS) -DJUDYNEXT -c $<

Judy1Prev.o: Judy1Prev.c
	@echo "        CC $<"
	$(CC) $(CFLAGS) -DJUDYPREV -c $<

Judy1PrevEmpty.o: Judy1PrevEmpty.c
	@echo "        CC $<"
	$(CC) $(CFLAGS) -DJUDYPREV -c $<

Judy1ByCount.o: Judy1ByCount.c
	@echo "        CC $<"
	$(CC) $(CFLAGS) -DNOSMARTJBB -DNOSMARTJBU -DNOSMARTJLB -c $<

Judy1Tables.o: Judy1Tables.c
	@echo "        CC $<"
	$(CC) $(CFLAGS) -c $<

Judy1Tables.c:
	$(CC) $(CFLAGS) -DJUDY1 Judy1TablesGen.c -o Judy1TablesGen
	./Judy1TablesGen

clean:
	rm -f *.o 
	rm -f Judy1Tables.c

