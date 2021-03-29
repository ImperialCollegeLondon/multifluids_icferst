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
FCFLAGS = -I../include -I/usr/include  -I/usr/include/python3.6m -DHAVE_NUMPY -I/home/aniketj/.local/lib/python3.6/site-packages/numpy/core/include -I/usr/lib/petscdir/3.8.3/include -I/usr/lib/petscdir/3.8.3/linux-gnu-c-opt/include -I/usr/include/hdf5/openmpi -I/home/aniketj/Documents/latest_icferst/include/ -DHAVE_PETSC  -I/usr/include/vtk-7.1 -I/usr/include/freetype2 -I/usr/include -I/usr/include/x86_64-linux-gnu -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/opal/mca/event/libevent2022/libevent -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/opal/mca/event/libevent2022/libevent/include -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/include/python3.6m -I/usr/include/hdf5/openmpi -I/usr/include/jsoncpp -I/usr/include/libxml2 -I/usr/include/tcl   -DHAVE_VTK -ffast-math -frecord-marker=4  -DHAVE_NUMPY -I/home/aniketj/.local/lib/python3.6/site-packages/numpy/core/include -fno-realloc-lhs -ffree-line-length-none -ffixed-line-length-none  -O3 -I/usr/lib/petscdir/3.8.3/include -I/usr/lib/petscdir/3.8.3/linux-gnu-c-opt/include -I/usr/include/hdf5/openmpi -I/home/aniketj/Documents/latest_icferst/include/ -fdefault-real-8 -fdefault-double-8

CC  = mpicc
CFLAGS  =  -I/usr/include/python3.6m -DHAVE_NUMPY -I/home/aniketj/.local/lib/python3.6/site-packages/numpy/core/include  -O3 -I../include

CXX = mpicxx
CXXFLAGS=  -I/usr/include/python3.6m -DHAVE_NUMPY -I/home/aniketj/.local/lib/python3.6/site-packages/numpy/core/include -I/usr/lib/petscdir/3.8.3/include -I/usr/lib/petscdir/3.8.3/linux-gnu-c-opt/include -I/usr/include/hdf5/openmpi -I/home/aniketj/Documents/latest_icferst/include/ -DHAVE_PETSC  -I/usr/include/vtk-7.1 -I/usr/include/freetype2 -I/usr/include -I/usr/include/x86_64-linux-gnu -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/opal/mca/event/libevent2022/libevent -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/opal/mca/event/libevent2022/libevent/include -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/include/python3.6m -I/usr/include/hdf5/openmpi -I/usr/include/jsoncpp -I/usr/include/libxml2 -I/usr/include/tcl   -DHAVE_VTK  -I/usr/include/python3.6m -DHAVE_NUMPY -I/home/aniketj/.local/lib/python3.6/site-packages/numpy/core/include  -O3 -I../include

LINKER = mpicxx  -Xlinker -export-dynamic -Wl,-O1 -Wl,-Bsymbolic-functions  

FLLINKER = mpif90  -Xlinker -export-dynamic -Wl,-O1 -Wl,-Bsymbolic-functions   $(FCFLAGS)

LIBS = -L../lib ./lib/libadaptivity.a /home/aniketj/Documents/latest_icferst/lib/libspatialindex.a /home/aniketj/Documents/latest_icferst/lib/libspud.a ./lib/libvtkfortran.a ./lib/libjudy.a -lmba2d -lzoltan -Wl,-rpath,/usr/lib/petscdir/3.8.3/linux-gnu-c-opt/lib -Wl,-rpath,/usr/lib/petscdir/3.8.3/linux-gnu-c-opt/lib -L/usr/lib/petscdir/3.8.3/linux-gnu-c-opt/lib -Wl,-rpath,/usr/lib/petscdir/3.8.3/linux-gnu-c-opt/lib -Wl,-rpath,/usr/lib/x86_64-linux-gnu/openmpi/lib -L/usr/lib/x86_64-linux-gnu/openmpi/lib -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/7 -L/usr/lib/gcc/x86_64-linux-gnu/7 -Wl,-rpath,/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -lpetsc -lexoIIv2c -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lnetcdf -lparmetis -lumfpack -lklu -lcholmod -lbtf -lccolamd -lcolamd -lcamd -lamd -lsuitesparseconfig -lml -lscalapack -lHYPRE -ltriangle -lhdf5_openmpi -lhdf5_openmpi_hl -lmetis -lchaco -lptesmumps -lptscotch -lptscotcherr -lesmumps -lscotch -lscotcherr -llapack -lblas -lctetgen -lX11 -lm -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lgfortran -lm -lgfortran -lm -lquadmath -lmpi_cxx -lstdc++ -lm -lrt -lm -lpthread -lz -Wl,-rpath,/usr/lib/x86_64-linux-gnu/openmpi/lib -L/usr/lib/x86_64-linux-gnu/openmpi/lib -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/7 -L/usr/lib/gcc/x86_64-linux-gnu/7 -Wl,-rpath,/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -Wl,-rpath,/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -ldl -lmpi -lgcc_s -lpthread -ldl -lnetcdff -lnetcdf  -llapack -ludunits2 -lpthread -lm -lstdc++  -L/usr/lib -lpython3.6m -Xlinker -export-dynamic -Wl,-O1 -Wl,-Bsymbolic-functions -lexpat                     -L/usr/lib -lz  -lexpat -lpthread -ldl  -lutil -lm  -L/usr/lib/gcc/x86_64-linux-gnu/7 -L/usr/lib/gcc/x86_64-linux-gnu/7/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/7/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/7/../../.. -lgfortran -lm -lquadmath  -L/usr/lib/gcc/x86_64-linux-gnu/7 -L/usr/lib/gcc/x86_64-linux-gnu/7/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/7/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/7/../../.. -lgfortran -lm -lquadmath -L./lib -Wl,-rpath,/usr/lib/x86_64-linux-gnu/hdf5/openmpi:/usr/lib/x86_64-linux-gnu/openmpi/lib /usr/lib/x86_64-linux-gnu/libfreetype.so /usr/lib/x86_64-linux-gnu/libz.so /usr/lib/x86_64-linux-gnu/libexpat.so /usr/lib/x86_64-linux-gnu/libvtkDomainsChemistryOpenGL2-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libjpeg.so /usr/lib/x86_64-linux-gnu/libpng.so /usr/lib/x86_64-linux-gnu/libtiff.so /usr/lib/x86_64-linux-gnu/libvtkFiltersGeneric-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkFiltersHyperTree-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkFiltersParallelDIY2-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkFiltersParallelFlowPaths-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkFiltersParallelGeometry-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkFiltersParallelImaging-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkFiltersParallelMPI-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkFiltersParallelStatistics-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkFiltersPoints-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkFiltersProgrammable-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkFiltersPython-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libpython3.6m.so /usr/lib/libvtkWrappingTools-7.1.a /usr/lib/x86_64-linux-gnu/libvtkFiltersReebGraph-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkFiltersSMP-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkFiltersSelection-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkFiltersTexture-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkFiltersVerdict-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkverdict-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libproj.so /usr/lib/x86_64-linux-gnu/libvtkIOAMR-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/hdf5/openmpi/libhdf5.so /usr/lib/x86_64-linux-gnu/libsz.so /usr/lib/x86_64-linux-gnu/libdl.so -lm /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_cxx.so /usr/lib/x86_64-linux-gnu/hdf5/openmpi/libhdf5_hl.so /usr/lib/x86_64-linux-gnu/libvtkIOEnSight-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libnetcdf_c++.so /usr/lib/x86_64-linux-gnu/libnetcdf.so /usr/lib/x86_64-linux-gnu/libvtkIOExport-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkRenderingGL2PSOpenGL2-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libgl2ps.so /usr/lib/x86_64-linux-gnu/libvtkIOFFMPEG-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkIOMovie-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libtheoraenc.so /usr/lib/x86_64-linux-gnu/libtheoradec.so /usr/lib/x86_64-linux-gnu/libogg.so /usr/lib/x86_64-linux-gnu/libvtkIOGDAL-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkIOGeoJSON-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libjsoncpp.so /usr/lib/x86_64-linux-gnu/libvtkIOImport-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkIOInfovis-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libxml2.so /usr/lib/x86_64-linux-gnu/libvtkIOMINC-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkIOMPIImage-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkIOMPIParallel-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkIOParallel-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkIONetCDF-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkIOMySQL-7.1.so.7.1.1 -lsqlite3 /usr/lib/x86_64-linux-gnu/libvtkIOODBC-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkIOPLY-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkIOParallelExodus-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkIOParallelLSDyna-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkIOParallelNetCDF-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkIOParallelXML-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkIOPostgreSQL-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkIOTecplotTable-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkIOVPIC-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkVPIC-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkIOVideo-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkIOXdmf2-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkxdmf2-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkImagingMorphological-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkImagingStatistics-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkImagingStencil-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkInfovisBoostGraphAlgorithms-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkInteractionImage-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkLocalExample-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkParallelMPI4Py-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkRenderingContextOpenGL2-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkRenderingExternal-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkRenderingFreeTypeFontConfig-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkRenderingImage-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkRenderingLOD-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkRenderingMatplotlib-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkRenderingParallel-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkRenderingParallelLIC-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkRenderingSceneGraph-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkRenderingVolumeAMR-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkRenderingVolumeOpenGL2-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkTestingGenericBridge-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkTestingIOSQL-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkTestingRendering-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkViewsContext2D-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkViewsGeovis-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkViewsInfovis-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkWrappingJava-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkDomainsChemistry-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkFiltersFlowPaths-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkIOGeometry-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkIOExodus-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkexoIIc-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libnetcdf_c++.so /usr/lib/x86_64-linux-gnu/libnetcdf.so /usr/lib/x86_64-linux-gnu/libvtkIOLSDyna-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkIOSQL-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/hdf5/openmpi/libhdf5.so /usr/lib/x86_64-linux-gnu/libsz.so /usr/lib/x86_64-linux-gnu/libdl.so -lm /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_cxx.so /usr/lib/x86_64-linux-gnu/hdf5/openmpi/libhdf5_hl.so /usr/lib/x86_64-linux-gnu/libxml2.so /usr/lib/x86_64-linux-gnu/libvtkWrappingPython36Core-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkPythonInterpreter-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libpython3.6m.so /usr/lib/x86_64-linux-gnu/libvtkFiltersParallel-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkParallelMPI-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkRenderingLICOpenGL2-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkFiltersAMR-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkParallelCore-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkIOLegacy-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkRenderingOpenGL2-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libGLEW.so /usr/lib/x86_64-linux-gnu/libSM.so /usr/lib/x86_64-linux-gnu/libICE.so /usr/lib/x86_64-linux-gnu/libX11.so /usr/lib/x86_64-linux-gnu/libXext.so /usr/lib/x86_64-linux-gnu/libXt.so /usr/lib/x86_64-linux-gnu/libvtkImagingMath-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkGeovisCore-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libproj.so /usr/lib/x86_64-linux-gnu/libvtkChartsCore-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkRenderingContext2D-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkFiltersImaging-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkInfovisLayout-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkInfovisCore-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkViewsCore-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkInteractionWidgets-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkFiltersHybrid-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkImagingGeneral-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkImagingSources-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkFiltersModeling-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkImagingHybrid-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkIOImage-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkDICOMParser-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkmetaio-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libz.so /usr/lib/x86_64-linux-gnu/libvtkInteractionStyle-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkFiltersExtraction-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkFiltersStatistics-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkImagingFourier-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkalglib-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkRenderingAnnotation-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkImagingColor-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkRenderingVolume-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkImagingCore-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkIOXML-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkIOXMLParser-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkIOCore-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkRenderingLabel-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkRenderingFreeType-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkRenderingCore-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkCommonColor-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkFiltersGeometry-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkFiltersSources-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkFiltersGeneral-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkCommonComputationalGeometry-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkFiltersCore-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkCommonExecutionModel-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkCommonDataModel-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkCommonTransforms-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkCommonMisc-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkCommonMath-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtkCommonSystem-7.1.so.7.1.1 /usr/lib/x86_64-linux-gnu/libvtksys-7.1.so.7.1.1 -ldl /usr/lib/x86_64-linux-gnu/libfreetype.so /usr/lib/x86_64-linux-gnu/libvtkCommonCore-7.1.so.7.1.1   -lblas  -L/usr/lib/gcc/x86_64-linux-gnu/7 -L/usr/lib/gcc/x86_64-linux-gnu/7/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/7/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/7/../../.. -lgfortran -lm -lquadmath 
FLIBS = -L../lib  -L/usr/lib/gcc/x86_64-linux-gnu/7 -L/usr/lib/gcc/x86_64-linux-gnu/7/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/7/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/7/../../.. -lgfortran -lm -lquadmath

# These objects are needed in each program as they provide the code related to
# diagnostic output and FLabort:
DEBUG_OBJECTS=../debug/Debug.o ../debug/Debug_Parameters.o

VERTICAL_INTEGRATION = ../bin/vertical_integration
VTKDIAGNOSTIC = ../bin/vtkdiagnostic
FLDIAGNOSTICS = ../bin/fldiagnostics
STREAMFUNCTION_2D = ../bin/streamfunction_2d
FLREDECOMP = ../bin/flredecomp
MESHCONV = ../bin/meshconv
TEST_LAPLACIAN= ../bin/test_laplacian
VISUALISE_ELEMENTS= ../bin/visualise_elements
GMSH2VTU= ../bin/gmsh2vtu
TRIANGLE2VTU= ../bin/triangle2vtu
VTU2GMSH= ../bin/vtu2gmsh
VTU2TRIANGLE= ../bin/vtu2triangle
PETSC_READNSOLVE = ../bin/petsc_readnsolve
TEST_PRESSURE_SOLVE = ../bin/test_pressure_solve
PROJECT_TO_CONTINUOUS = ../bin/project_to_continuous
UNIFIEDMESH = ../bin/unifiedmesh
LIBFLUIDITY=-l$(FLUIDITY)
LIBFEMTOOLS=-lfemtools
LIBADAPT=-ladaptivity
CHECKMESH=../bin/checkmesh
PROJECT_VTU=../bin/project_vtu
PROBE_VTU=../bin/probe_vtu
SUPERMESH_DIFFERENCE=../bin/supermesh_difference
DIFFERENTIATE_VTU=../bin/differentiate_vtu
VTU_BINS=../bin/vtu_bins
FLADAPT=../bin/fladapt
VTKPROJECTION=../bin/vtk_projection
PERIODISE=../bin/periodise

BINARIES = $(VTKDIAGNOSTIC)		\
  $(FLDIAGNOSTICS) $(FLREDECOMP) $(PETSC_READNSOLVE)			\
  $(VERTICAL_INTEGRATION) $(EVMERGE) $(FLADAPT) $(CHECKMESH)		\
  $(PROJECT_VTU) $(PROBE_VTU) $(SUPERMESH_DIFFERENCE)			\
  $(DIFFERENTIATE_VTU) $(VTU_BINS) $(GMSH2VTU)	        		\
  $(PROJECT_TO_CONTINUOUS) $(MESHCONV) $(VTU2GMSH)		        \
  $(TEST_PRESSURE_SOLVE) $(UNIFIEDMESH) $(VTKPROJECTION) $(PERIODISE)	\
  $(STREAMFUNCTION_2D)							\
  $(TEST_LAPLACIAN) 

# Work around Intel compiler bug in visualise elements by simply not
# building it.
ifeq (yes,yes)
BINARIES += $(VISUALISE_ELEMENTS)
endif


.SUFFIXES: .f90 .F90 .c .cpp .o .a 

%.o:	%.f90
	$(FC) $(FCFLAGS) -c $<
%.o:	%.F90
	$(FC) $(FCFLAGS) -c $<
%.o:	%.c
	$(CC) $(CFLAGS) -c $<
%.o:	%.cpp
	$(CXX) $(CXXFLAGS) -c $<

default: $(BINARIES) scripts

scripts:
	mkdir -p ../bin
	for i in *.py ; do cp $$i ../bin/`basename $$i .py` ; done
	for i in *.sh ; do cp $$i ../bin/`basename $$i .sh` ; done

install-version-info: version-info
	@mkdir -p ../bin
	@chmod +x $<
	@cp $< ../bin/$<

.PHONY: install-version-info

.PHONY: install_python

install_python:
	python3 setup.py bdist --format=tar
	tar -x -C /$(DESTDIR) -f dist/*.tar

$(GMSH2VTU): gmsh2vtu.o gmsh2vtu_main.o lib/
	$(LINKER) -o $@ $(filter %.o,$^) -l$(FLUIDITY) $(LIBS)

$(TRIANGLE2VTU): triangle2vtu.o triangle2vtu_main.o lib/
	$(LINKER) -o $@ $(filter %.o,$^) -l$(FLUIDITY) $(LIBS)

$(VTU2GMSH): vtu2gmsh.o vtu2gmsh_main.o lib/
	$(LINKER) -o $@ $(filter %.o,$^) -l$(FLUIDITY) $(LIBS)

$(VTU2TRIANGLE): vtu2triangle.o vtu2triangle_main.o lib/
	$(LINKER) -o $@ $(filter %.o,$^) -l$(FLUIDITY) $(LIBS)

$(PROJECT_TO_CONTINUOUS): project_to_continuous.o project_to_continuous_main.o
	$(LINKER) -o $@ $(filter %.o,$^) -l$(FLUIDITY) $(LIBS)

$(TRIANGLE_TEST_MESH): triangle_test_mesh.o
	$(LINKER) -o $(TRIANGLE_TEST_MESH) $^ $(FLIBS) $(LIBFLUIDITY) $(LIBS)

$(TEST_LAPLACIAN): test_laplacian.o
	$(FLLINKER) -o $@ $(filter %.o,$^) -l$(FLUIDITY) $(LIBS)

$(VISUALISE_ELEMENTS): visualise_elements.o
	$(FLLINKER) $(CXXFLAGS) -o $@ $(filter %.o,$^) -l$(FLUIDITY) $(LIBS)

$(VERTICAL_INTEGRATION): Vertical_Integration.o Vertical_Integration_main.o lib/
	$(LINKER) $(CXXFLAGS) -o $@ $(filter %.o,$^) -l$(FLUIDITY) $(LIBS)

$(PETSC_READNSOLVE): petsc_readnsolve_main.o petsc_readnsolve.o \
  lib/
	$(LINKER) $(CXXFLAGS) -o $@ $(filter %.o,$^) -l$(FLUIDITY) $(LIBS)

$(TEST_PRESSURE_SOLVE): test_pressure_solve_main.o test_pressure_solve.o \
  lib/
	$(LINKER) $(CXXFLAGS) -o $@ $(filter %.o,$^) -l$(FLUIDITY) $(LIBS)

# we need a symlink to ../lib, as $(LIBS) contains relative references to
# libraries (i.e. ./lib/libadaptivity.a instead of -ladaptivity)
lib/:
	ln -sf ../lib

$(VTKDIAGNOSTIC): vtkdiagnostic.cpp lib/
	$(LINKER) $(CXXFLAGS) -I../include -o $(VTKDIAGNOSTIC) vtkdiagnostic.cpp $(LIBS)

$(VTKPROJECTION): vtkprojection.cpp lib/
	$(LINKER) $(CXXFLAGS) -I../include -o $(VTKPROJECTION) vtkprojection.cpp -L../lib/ -l$(FLUIDITY) $(LIBS)

$(PERIODISE): periodise.o lib/
	$(FLLINKER) -o $@ $(filter %.o,$^) -l$(FLUIDITY) $(LIBS)

$(FLDIAGNOSTICS): fldiagnostics.o fldiagnostics_main.o lib/
	$(LINKER) -o $@ $(filter %.o,$^) -l$(FLUIDITY) $(LIBS)

$(STREAMFUNCTION_2D): Streamfunction_2D.o Streamfunction_2D_main.o lib/
	$(LINKER) -o $@ $(filter %.o,$^) -l$(FLUIDITY) $(LIBS)

$(FLREDECOMP): Flredecomp.o Flredecomp_main.o lib/
	$(LINKER) -o $@ $(filter %.o,$^) -l$(FLUIDITY) $(LIBS)

$(MESHCONV): Meshconv.o Meshconv_main.o lib/
	$(LINKER) -o $@ $(filter %.o,$^) -l$(FLUIDITY) $(LIBS)

$(CHECKMESH): Checkmesh.o Checkmesh_main.o lib/
	$(LINKER) -o $@ $(filter %.o,$^) -l$(FLUIDITY) $(LIBS)

$(PROJECT_VTU): Project_Vtu.o Project_Vtu_main.o lib/
	$(LINKER) -o $@ $(filter %.o,$^) -l$(FLUIDITY) $(LIBS)

$(PROBE_VTU): Probe_Vtu.o Probe_Vtu_main.o lib/
	$(LINKER) -o $@ $(filter %.o,$^) -l$(FLUIDITY) $(LIBS)

$(SUPERMESH_DIFFERENCE): Supermesh_Difference.o Supermesh_Difference_main.o lib/
	$(LINKER) -o $@ $(filter %.o,$^) -l$(FLUIDITY) $(LIBS)

$(DIFFERENTIATE_VTU): Differentiate_Vtu.o Differentiate_Vtu_main.o lib/
	$(LINKER) -o $@ $(filter %.o,$^) -l$(FLUIDITY) $(LIBS)

$(VTU_BINS): Vtu_Bins.o Vtu_Bins_main.o lib/
	$(LINKER) -o $@ $(filter %.o,$^) -l$(FLUIDITY) $(LIBS)

$(FLADAPT): Fladapt.o Fladapt_main.o lib/
	$(LINKER) -o $@ $(filter %.o,$^) -l$(FLUIDITY) $(LIBS)

$(UNIFIEDMESH): unifiedmesh.o unifiedmesh_main.o lib/
	$(LINKER) -o $@ $(filter %.o,$^) -l$(FLUIDITY) $(LIBS)

clean:
	rm -f *.o *.mod
	rm -rf build dist fluidity.egg-info 
	rm -f lib

distclean:
	rm -f version-info

