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
FCFLAGS = -I../../include -ffast-math -frecord-marker=4  -DHAVE_NUMPY -I/usr/lib/python3/dist-packages/numpy/core/include -fno-realloc-lhs -ffree-line-length-none -ffixed-line-length-none  -O3 -I/media/hh210/data2/Petsc/petsc-3.19.5/include -I/media/hh210/data2/Petsc/petsc-3.19.5/arch-linux-c-debug/include -I/usr/include/python3.12 -fdefault-real-8 -fdefault-double-8 -I../

CC  = mpicc
CFLAGS  =  -I/usr/include/python3.12 -I/usr/include/python3.12 -DHAVE_NUMPY -I/usr/lib/python3/dist-packages/numpy/core/include -DFLUIDITY_PYTHON_MODULE_PATH=\"/media/hh210/data2/ICFEMRES/ICFERST_petsc3.19/python\"  -O3 -I/media/hh210/data2/Petsc/petsc-3.19.5/include -I/media/hh210/data2/Petsc/petsc-3.19.5/arch-linux-c-debug/include -I../../include

CXX = mpicxx -I../../include
CXXFLAGS=  -I/usr/include/python3.12 -I/usr/include/python3.12 -DHAVE_NUMPY -I/usr/lib/python3/dist-packages/numpy/core/include -I/media/hh210/data2/Petsc/petsc-3.19.5/include -I/media/hh210/data2/Petsc/petsc-3.19.5/arch-linux-c-debug/include -DHAVE_PETSC  -I/usr/include/vtk-9.1 -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -std=c++11 -I../include -I../../include 

LDFLAGS  =   

LIBS = ./lib/libadaptivity.a ./lib/libvtkfortran.a /media/hh210/data2/ICFEMRES/ICFERST_petsc3.19/lib/libspatialindex.a /media/hh210/data2/ICFEMRES/ICFERST_petsc3.19/lib/libspud.a ./lib/libjudy.a -lmba2d -lparmetis -lmetis -L/media/hh210/data2/Petsc/petsc-3.19.5/arch-linux-c-debug/lib -lpetsc -lnetcdff -lnetcdf  -llapack -ludunits2 -lpthread -lm -lstdc++  -L/usr/lib/python3.12/config-3.12-x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -lpython3.12 -ldl  -lm  -ldl -lpthread -lutil -lm  -L/usr/lib/gcc/x86_64-linux-gnu/13 -L/usr/lib/gcc/x86_64-linux-gnu/13/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/13/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/13/../../.. -lgfortran -lm -lquadmath  -L/usr/lib/gcc/x86_64-linux-gnu/13 -L/usr/lib/gcc/x86_64-linux-gnu/13/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/13/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/13/../../.. -lgfortran -lm -lquadmath -L./lib  -lvtkCommonCore-9.1 -lvtkCommonDataModel-9.1 -lvtkCommonExecutionModel-9.1 -lvtkFiltersCore-9.1 -lvtkFiltersVerdict-9.1 -lvtkIOCore-9.1 -lvtkIOXML-9.1 -lvtkIOGeometry-9.1 -lvtkIOLegacy-9.1 -lvtkIOParallelXML-9.1 -lvtkParallelMPI-9.1 -lvtkIOXMLParser-9.1 -lvtkCommonSystem-9.1 -lvtkCommonMisc-9.1 -lvtkCommonMath-9.1 -lvtksys-9.1 -ldl -lpthread -L/usr/lib/x86_64-linux-gnu/openmpi/lib -lmpi_cxx -lmpi 
FLIBS =  -L/usr/lib/gcc/x86_64-linux-gnu/13 -L/usr/lib/gcc/x86_64-linux-gnu/13/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/13/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/13/../../.. -lgfortran -lm -lquadmath
LIBFLUIDITY=../../lib/lib$(FLUIDITY).a

# the test binaries NOT to be built
DISABLED_TESTS=test_anisotropic_adaptivity test_anisotropic_adaptivity_two test_gradation test_isotropic_gradation \
  test_annulus_gradation \
  test_directional_gradation test_directional_gradation_annulus test_pseudo2d_gradation \
  test_assemble_metric test_anisotropically_bounded_metric test_backstep \
  test_quadratic_interpolation \
  test_recovery_estimator test_elementwise_error_adaptivity \
  test_boundary_layer_adaptivity test_boundary_layer_adaptivity_two \
  test_anisotropic_bounds_equivalence \
  compute_interpolation_error_noadapt compute_interpolation_error_adapt \
  compute_enstrophy_goal compute_temperature_goal compute_hessian_error \
  compute_chimney_adapt compute_les_temp_goal compute_les_velocity_goal \
  compute_mesh_conformity compute_driven_cavity_adapt test_mba compute_residual \
  test_metric_advection test_geometric_constraints  \
  compute_compare_interpolation compute_compare_interpolation_3d \
  test_anisotropic_gradation compute_anisotropic_gradation \
  test_anisotropic_zz_3d test_anisotropic_zz


# the test programs to be built:
TEST_BINARIES=$(filter-out $(DISABLED_TESTS), $(basename $(wildcard *.F90)))

unittest: $(TEST_BINARIES)
	ln -sf $(addprefix ../../error_measures/tests/,$(TEST_BINARIES)) ../../bin/tests/

.SUFFIXES: .f90 .F90 .c .cpp .o .a $(.SUFFIXES)

%.o:	%.f90
	$(FC) $(FCFLAGS) -c $<
%.o:	%.F90
	$(FC) $(FCFLAGS) -c $<
%.o:	%.c
	$(CC) $(CFLAGS) -c $<
%.o:	%.cpp
	$(CXX) $(CXXFLAGS) -c $<

# creates a TESTNAME_main.o from test_main.cpp which calls the subroutine
# TESTNAME, that should be a subroutine in TESTNAME.F90
%_main.o:
	$(CXX) $(CXXFLAGS) -D TESTNAME=$(subst _main.o,,$@)_ -o $@ -c test_main.cpp 

# link this TESTNAME_main.o with TESTNAME.o from TESTNAME.F90
%: %_main.o %.o lib/
	$(CXX) $(LDFLAGS) -o $@ $(filter %.o,$^) $(LIBFLUIDITY) $(LIBS) 

# make a temp sym-link to the directory containing fluidity and adaptivity
lib/:
	ln -s ../../lib

clean:
	rm -f $(TEST_BINARIES)
	rm -f *.o *.mod
	rm -f lib
