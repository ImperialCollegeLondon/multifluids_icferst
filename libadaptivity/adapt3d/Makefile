# Copyright (C) 2006 Imperial College London and others.
#
# Please see the AUTHORS file in the main source directory for a full list
# of copyright holders.
#
# Adrian Umpleby
# Applied Modelling and Computation Group
# Department of Earth Science and Engineering
# Imperial College London
#
# adrian@Imperial.ac.uk
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
# USA

SHELL = /bin/bash

LIB 	= ../lib/libadaptivity.a

CPPFLAGS =  -I/usr/include/vtk-6.3 -I/usr/include/freetype2 -I/usr/include -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/opal/mca/event/libevent2022/libevent -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/opal/mca/event/libevent2022/libevent/include -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/include/python2.7 -I/usr/include/x86_64-linux-gnu -I/usr/include/hdf5/openmpi -I/usr/include/libxml2 -I/usr/include/jsoncpp -I/usr/include/tcl   -DHAVE_VTK -I./include -I../include

FC	= mpif90
FCFLAGS	= $(CPPFLAGS) -ffast-math -frecord-marker=4    -ffree-line-length-none -ffixed-line-length-none  -O3 -fdefault-real-8

CXX	= mpicxx
CXXFLAGS= $(CPPFLAGS)     -O3

LINKER	= mpicxx
LFLAGS	=  -llapack -lblas -lpthread -lm -lstdc++  -Wl,-rpath,/usr/lib/x86_64-linux-gnu/hdf5/openmpi:/usr/lib/x86_64-linux-gnu/openmpi/lib /usr/lib/x86_64-linux-gnu/libfreetype.so /usr/lib/x86_64-linux-gnu/libz.so /usr/lib/x86_64-linux-gnu/libvtkDomainsChemistry-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libexpat.so /usr/lib/x86_64-linux-gnu/libvtkFiltersGeneric-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkFiltersHyperTree-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkFiltersParallelFlowPaths-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkFiltersParallelGeometry-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkFiltersParallelImaging-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkFiltersParallelMPI-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkFiltersParallelStatistics-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkFiltersProgrammable-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkFiltersPython-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libpython2.7.so /usr/lib/libvtkWrappingTools-6.3.a /usr/lib/x86_64-linux-gnu/libvtkFiltersReebGraph-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkFiltersSMP-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkFiltersSelection-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkFiltersTexture-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkFiltersVerdict-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkverdict-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libjpeg.so /usr/lib/x86_64-linux-gnu/libpng.so /usr/lib/x86_64-linux-gnu/libtiff.so /usr/lib/x86_64-linux-gnu/libproj.so /usr/lib/x86_64-linux-gnu/libvtkIOAMR-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/hdf5/openmpi/libhdf5.so /usr/lib/x86_64-linux-gnu/libsz.so /usr/lib/x86_64-linux-gnu/libdl.so -lm /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so /usr/lib/x86_64-linux-gnu/libvtkIOEnSight-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libnetcdf_c++.so /usr/lib/x86_64-linux-gnu/libnetcdf.so /usr/lib/x86_64-linux-gnu/libvtkIOExport-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkRenderingGL2PS-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkRenderingContextOpenGL-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libgl2ps.so /usr/lib/x86_64-linux-gnu/libvtkIOFFMPEG-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkIOMovie-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libtheoraenc.so /usr/lib/x86_64-linux-gnu/libtheoradec.so /usr/lib/x86_64-linux-gnu/libogg.so /usr/lib/x86_64-linux-gnu/libvtkIOGDAL-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkIOGeoJSON-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkIOImport-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkIOInfovis-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libxml2.so /usr/lib/x86_64-linux-gnu/libvtkIOMINC-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkIOMPIImage-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkIOMPIParallel-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkIOParallel-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkIONetCDF-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libjsoncpp.so /usr/lib/x86_64-linux-gnu/libvtkIOMySQL-6.3.so.6.3.0 -lsqlite3 /usr/lib/x86_64-linux-gnu/libvtkIOODBC-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkIOPLY-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkIOParallelExodus-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkIOParallelLSDyna-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkIOParallelNetCDF-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkIOParallelXML-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkIOPostgreSQL-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkIOVPIC-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkVPIC-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkIOVideo-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkIOXdmf2-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkxdmf2-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkImagingMath-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkImagingMorphological-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkImagingStatistics-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkImagingStencil-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkInteractionImage-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkLocalExample-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkParallelMPI4Py-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkRenderingExternal-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkRenderingFreeTypeFontConfig-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkRenderingImage-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkRenderingLOD-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkRenderingMatplotlib-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkRenderingParallel-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkRenderingParallelLIC-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkRenderingVolumeAMR-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkRenderingVolumeOpenGL-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkTestingGenericBridge-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkTestingIOSQL-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkTestingRendering-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkViewsContext2D-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkViewsGeovis-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkViewsInfovis-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkWrappingJava-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkFiltersFlowPaths-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkIOExodus-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkexoIIc-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libnetcdf_c++.so /usr/lib/x86_64-linux-gnu/libnetcdf.so /usr/lib/x86_64-linux-gnu/libvtkIOLSDyna-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/hdf5/openmpi/libhdf5.so /usr/lib/x86_64-linux-gnu/libsz.so /usr/lib/x86_64-linux-gnu/libdl.so -lm /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so /usr/lib/x86_64-linux-gnu/libxml2.so /usr/lib/x86_64-linux-gnu/libvtkWrappingPython27Core-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkPythonInterpreter-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libpython2.7.so /usr/lib/x86_64-linux-gnu/libvtkFiltersParallel-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkParallelMPI-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkRenderingLIC-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkFiltersAMR-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkParallelCore-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkIOLegacy-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkRenderingOpenGL-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libGLU.so /usr/lib/x86_64-linux-gnu/libSM.so /usr/lib/x86_64-linux-gnu/libICE.so /usr/lib/x86_64-linux-gnu/libX11.so /usr/lib/x86_64-linux-gnu/libXext.so /usr/lib/x86_64-linux-gnu/libXt.so /usr/lib/x86_64-linux-gnu/libvtkIOSQL-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkGeovisCore-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkIOXML-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkIOGeometry-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkIOXMLParser-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libproj.so /usr/lib/x86_64-linux-gnu/libvtkChartsCore-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkRenderingContext2D-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkFiltersImaging-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkInfovisLayout-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkInfovisBoostGraphAlgorithms-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkInfovisCore-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkViewsCore-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkInteractionWidgets-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkFiltersHybrid-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkImagingGeneral-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkImagingSources-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkFiltersModeling-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkImagingHybrid-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkIOImage-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkDICOMParser-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkIOCore-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkmetaio-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libz.so /usr/lib/x86_64-linux-gnu/libvtkInteractionStyle-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkRenderingAnnotation-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkImagingColor-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkRenderingVolume-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkRenderingLabel-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkRenderingFreeType-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkRenderingCore-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkCommonColor-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkFiltersExtraction-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkFiltersStatistics-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkImagingFourier-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkImagingCore-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkalglib-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkFiltersGeometry-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkFiltersSources-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkFiltersGeneral-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkFiltersCore-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkCommonExecutionModel-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkCommonComputationalGeometry-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkCommonDataModel-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkCommonMisc-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkCommonTransforms-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkCommonMath-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtkCommonSystem-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libvtksys-6.3.so.6.3.0 -ldl /usr/lib/x86_64-linux-gnu/libvtkftgl-6.3.so.6.3.0 /usr/lib/x86_64-linux-gnu/libfreetype.so /usr/lib/x86_64-linux-gnu/libGL.so /usr/lib/x86_64-linux-gnu/libvtkCommonCore-6.3.so.6.3.0   -L/usr/lib/gcc/x86_64-linux-gnu/7 -L/usr/lib/gcc/x86_64-linux-gnu/7/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/7/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/7/../../.. -lgfortran -lm -lquadmath -L/usr/lib/gcc/x86_64-linux-gnu/7 -L/usr/lib/gcc/x86_64-linux-gnu/7/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/7/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/7/../../.. -lgfortran -lm -lquadmath -L./lib -L/usr/lib/gcc/x86_64-linux-gnu/7 -L/usr/lib/gcc/x86_64-linux-gnu/7/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/7/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/7/../../.. -lgfortran -lm -lquadmath

MAKE	= make
INSTALL = /usr/bin/install -c

VPATH	= ./src

OBJS = Debug.o AdaptProgress.o write_log.o adapt1.o  addedg.o  addele.o  addfac.o  addnd.o   adfrnd.o \
		adptvy.o  angedg.o  angfac.o  chkfac.o  chkgmy.o  cmlced.o \
		cmneds.o  cmnnod.o  cnctel.o  cnctfc.o  dofacs.o  edg2fc.o \
		edg3nd.o  edg4nd.o  edgdst.o  edgeig.o  edgels.o  edgfnc.o \
		edginf.o  edgtst.o  elmdiv.o  elmeig.o  elmfnc.o  elminf.o \
		elmtst.o  elnods.o  elstst.o  facinf.o  facnrm.o  fc2edg.o \
		fchsnd.o  fixeig.o  fndedg.o  fndelm.o  fndsrf.o  getflg.o \
		getmsh.o  getsiz.o  gtfrfl.o  gtndfl.o  hasedg.o  adaptmem.o \
		hasnod.o  inipwr.o  inirls.o  initim.o  initls.o  intfld.o \
		intrpl.o  lenedg.o  lisbig.o  lstels.o  lstfac.o  lstnds.o \
		lstout.o  minang.o  mk1ele.o  mkadpt.o  mkecls.o  mkfixd.o \
		mkfrnt.o  mkncls.o  mknels.o  mshout.o  ndcnct.o  nodeig.o \
		nodfcs.o  nodmov.o  nodtst.o  nxtid.o   oppedg.o  ordnds.o \
		qlyout.o  rmfrnd.o  rmvedg.o  rmvele.o  rmvfac.o  rmvnd.o \
		rmvnxt.o  rspnse.o  second.o  sectlp.o  setflg.o  setgmy.o \
		shredg.o  shwhst.o  shwtim.o  spledg.o  stchfr.o  stfrfl.o \
		stfrgm.o  stndfl.o  mtetin.o  tetvol.o  undstt.o  vals3d.o \
		invrse.o  chkint.o  assval.o  nwcnel.o  wchfac.o  chgeds.o \
		Adaptivity.o expected_elements.o Flag_Handling.o

.SUFFIXES: .F .F90 .c .o .a

.F.o:
	$(FC) $(FCFLAGS) -c $<

.F90.o:
	$(FC) $(FCFLAGS) -c $<

.cpp.a:
	$(CXX) $(CXXFLAGS) -c $<

$(LIB): $(OBJS)
	@mkdir -p ../lib
	@$(AR) $(ARFLAGS) $(LIB) $(OBJS)

include/fdebug_adapt.h: include/confdefs.h include/ewrite.h

Debug.o: include/fdebug_adapt.h write_log.o

adapt1.o :	adapt1.F AdaptProgress.o write_log.o include/blknew.i include/blkbig.i include/blkorg.i include/blknxt.i include/blkerr.i include/blktim.i

addedg.o :	addedg.F include/blknew.i include/blkbig.i include/blkerr.i

addele.o :	addele.F write_log.o include/blknew.i include/blkbig.i include/blkerr.i include/blktim.i

addfac.o :	addfac.F include/blkfrt.i include/blknwf.i include/blkerr.i

addnd.o :	addnd.F include/blknew.i include/blkbig.i include/blkorg.i include/blkerr.i

adfrnd.o :	adfrnd.F include/blknew.i include/blkbig.i include/blkerr.i

adptvy.o :	adptvy.F AdaptProgress.o write_log.o include/blkprm.i include/blknew.i include/blkerr.i

angedg.o :	angedg.F include/blknew.i include/blkbig.i include/blkerr.i

angfac.o :	angfac.F include/blknew.i include/blkbig.i include/blkerr.i

chkfac.o :	chkfac.F include/blknew.i include/blkbig.i include/blkfrt.i include/blkerr.i

chkgmy.o :	chkgmy.F write_log.o include/blknew.i include/blkbig.i include/blkerr.i

cmlced.o :	cmlced.F include/blknew.i include/blkbig.i include/blkerr.i

cmneds.o :	cmneds.F include/blknew.i include/blkbig.i include/blkerr.i

cmnnod.o :	cmnnod.F include/blknew.i include/blkbig.i include/blkerr.i

cnctel.o :	cnctel.F include/blknew.i include/blkbig.i include/blkerr.i

nwcnel.o :	nwcnel.F include/blknew.i include/blkbig.i include/blkerr.i

wchfac.o :	wchfac.F include/blknew.i include/blkbig.i include/blkerr.i

chgeds.o :	chgeds.F include/blknew.i include/blkbig.i include/blkerr.i

cnctfc.o :	cnctfc.F include/blkfrt.i include/blkerr.i

dofacs.o :	dofacs.F include/blknew.i include/blkbig.i include/blkfrt.i include/blkerr.i

edg2fc.o :	edg2fc.F include/blknew.i include/blkbig.i include/blkerr.i include/blktim.i

edg3nd.o :	edg3nd.F include/blknew.i include/blkbig.i include/blkerr.i

edg4nd.o :	edg4nd.F include/blknew.i include/blkbig.i include/blkerr.i

edgdst.o :	edgdst.F include/blknew.i include/blkbig.i include/blkerr.i

edgeig.o :	edgeig.F include/blknew.i include/blkbig.i include/fixeig.i

edgels.o :	edgels.F include/blknew.i include/blkbig.i include/blkerr.i include/blktim.i

edgfnc.o :	edgfnc.F include/blknew.i include/blkbig.i include/blkerr.i

edginf.o :	edginf.F include/blknew.i include/blkbig.i include/blkerr.i

edgtst.o :	edgtst.F include/blknew.i include/blkbig.i include/blkorg.i include/blkerr.i include/blktim.i

elmdiv.o :	elmdiv.F include/blknew.i include/blkbig.i include/blkerr.i include/blkorg.i

elmeig.o :	elmeig.F include/blknew.i include/blkbig.i include/fixeig.i

elmfnc.o :	elmfnc.F write_log.o include/blknew.i include/blkbig.i include/blkerr.i

elminf.o :	elminf.F include/blknew.i include/blkbig.i include/blkerr.i

elmtst.o :	elmtst.F include/blknew.i include/blkbig.i include/blkerr.i include/blkorg.i

elnods.o :	elnods.F include/blknew.i include/blkbig.i include/blkerr.i

elstst.o :	elstst.F include/blknew.i include/blkbig.i include/blkerr.i include/blktim.i

facinf.o :	facinf.F include/blkfrt.i include/blknwf.i include/blkerr.i

facnrm.o :	facnrm.F include/blknew.i include/blkbig.i include/blkerr.i

fc2edg.o :	fc2edg.F write_log.o include/blknew.i include/blkbig.i include/blkerr.i include/blktim.i

fchsnd.o :	fchsnd.F include/blknew.i include/blkbig.i include/blkfrt.i include/blkerr.i

fixeig.o :	fixeig.F include/fixeig.i

fndedg.o :	fndedg.F

fndelm.o :	fndelm.F include/blkorg.i include/blkerr.i

fndsrf.o :	fndsrf.F

getflg.o :	getflg.F include/blknew.i include/blkbig.i include/blkerr.i include/blkpwr.i Flag_Handling.o

getmsh.o :	getmsh.F

getsiz.o :	getsiz.F include/blknew.i include/blkbig.i

gtfrfl.o :	gtfrfl.F include/blkfrt.i include/blknwf.i include/blkerr.i include/blkpwr.i

gtndfl.o :	gtndfl.F include/blknew.i include/blkbig.i include/blkerr.i include/blkpwr.i Flag_Handling.o

hasedg.o :	hasedg.F include/blknew.i include/blkbig.i include/blkerr.i

adaptmem.o	:	adaptmem.F

hasnod.o :	hasnod.F include/blknew.i include/blkbig.i include/blkerr.i

inipwr.o :	inipwr.F include/blkpwr.i

inirls.o :	inirls.F include/blkerr.i

initim.o :	initim.F include/blktim.i

initls.o :	initls.F include/blkerr.i

intfld.o :	intfld.F include/blknew.i include/blkbig.i include/blkerr.i

intrpl.o :	intrpl.F write_log.o include/blkerr.i

lenedg.o :	lenedg.F include/blknew.i include/blkbig.i include/blkerr.i

lisbig.o :	lisbig.F include/blknew.i include/blkbig.i include/blkerr.i

lstels.o :	lstels.F include/blknew.i include/blkbig.i include/blkerr.i

lstfac.o :	lstfac.F include/blkfrt.i include/blknwf.i include/blkerr.i

lstnds.o :	lstnds.F include/blknew.i include/blkbig.i include/blkerr.i

lstout.o :	lstout.F include/blknew.i include/blkbig.i include/blkerr.i

minang.o :	minang.F include/blknew.i include/blkbig.i include/blkfrt.i include/blkerr.i

mk1ele.o :	mk1ele.F include/blknew.i include/blkbig.i include/blkfrt.i include/blkerr.i

mkadpt.o :	mkadpt.F write_log.o include/blknew.i include/blkbig.i include/blkerr.i

mkecls.o :	mkecls.F

mkfixd.o :	mkfixd.F write_log.o include/blknew.i include/blkbig.i include/blkerr.i

mkfrnt.o :	mkfrnt.F include/blknew.i include/blkbig.i include/blkfrt.i include/blknwf.i include/blkerr.i

mkncls.o :	mkncls.F

mknels.o :	mknels.F

mshout.o :	mshout.F

ndcnct.o :	ndcnct.F include/blknew.i include/blkbig.i include/blkerr.i include/blktim.i

nodeig.o :	nodeig.F include/blknew.i include/blkbig.i include/fixeig.i

nodfcs.o :	nodfcs.F include/blknew.i include/blkbig.i include/blkfrt.i include/blkerr.i

nodmov.o :	nodmov.F write_log.o include/blknew.i include/blkbig.i include/blkorg.i include/blkerr.i include/blktim.i

nodtst.o :	nodtst.F include/blknew.i include/blkbig.i include/blkorg.i include/blkerr.i include/blktim.i

nxtid.o :	nxtid.F include/blkerr.i include/blknxt.i

oppedg.o :	oppedg.F include/blknew.i include/blkbig.i include/blkerr.i

ordnds.o :	ordnds.F include/blknew.i include/blkbig.i include/blkerr.i

qlyout.o :	qlyout.F include/blknew.i include/blkbig.i

rmfrnd.o :	rmfrnd.F include/blknew.i include/blkbig.i include/blkerr.i

rmvedg.o :	rmvedg.F include/blknew.i include/blkbig.i include/blknxt.i include/blkerr.i

rmvele.o :	rmvele.F write_log.o include/blknew.i include/blkbig.i include/blknxt.i include/blkerr.i include/blktim.i

rmvfac.o :	rmvfac.F include/blkfrt.i include/blknwf.i include/blkerr.i

rmvnd.o :	rmvnd.F include/blknew.i include/blkbig.i include/blkerr.i

rmvnxt.o :	rmvnxt.F include/blknxt.i

rspnse.o :	rspnse.F

second.o :	second.F

sectlp.o :	sectlp.F

setflg.o :	setflg.F include/blknew.i include/blkbig.i include/blkerr.i include/blkpwr.i Flag_Handling.o

setgmy.o :	setgmy.F write_log.o include/blknew.i include/blkbig.i include/blkerr.i

shredg.o :	shredg.F write_log.o include/blknew.i include/blkbig.i include/blkorg.i include/blkerr.i include/blktim.i

shwhst.o :	shwhst.F write_log.o include/blknew.i include/blkbig.i include/blkerr.i

shwtim.o :	shwtim.F write_log.o include/blktim.i

spledg.o :	spledg.F write_log.o include/blknew.i include/blkbig.i include/blkorg.i include/blkerr.i include/blktim.i

stchfr.o :	stchfr.F include/blknew.i include/blkbig.i include/blkfrt.i include/blknwf.i include/blkerr.i

stfrfl.o :	stfrfl.F include/blkfrt.i include/blknwf.i include/blkerr.i include/blkpwr.i

stfrgm.o :	stfrgm.F include/blknew.i include/blkbig.i include/blkfrt.i include/blkerr.i

stndfl.o :	stndfl.F include/blknew.i include/blkbig.i include/blkerr.i include/blkpwr.i Flag_Handling.o

mtetin.o :	mtetin.F

tetvol.o :	tetvol.F

undstt.o :	undstt.F

vals3d.o :	vals3d.F

invrse.o :	invrse.F

chkint.o :	chkint.F write_log.o include/blknew.i include/blkbig.i include/blkerr.i

assval.o :	assval.F

Flag_Handling.o: include/blkerr.i include/blkprm.i Debug.o

clean:
	rm -rf *.o *.mod

