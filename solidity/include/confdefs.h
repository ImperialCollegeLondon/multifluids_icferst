#define PACKAGE_NAME ""
#define PACKAGE_TARNAME ""
#define PACKAGE_VERSION ""
#define PACKAGE_STRING ""
#define PACKAGE_BUGREPORT ""
#define PACKAGE_URL ""
#define STDC_HEADERS 1
#define HAVE_SYS_TYPES_H 1
#define HAVE_SYS_STAT_H 1
#define HAVE_STDLIB_H 1
#define HAVE_STRING_H 1
#define HAVE_MEMORY_H 1
#define HAVE_STRINGS_H 1
#define HAVE_INTTYPES_H 1
#define HAVE_STDINT_H 1
#define HAVE_UNISTD_H 1
#define SIZEOF_LONG_INT 8
#define LONG_64_BITS 1
#define F77_FUNC(name,NAME) name ## _
#define F77_FUNC_(name,NAME) name ## _
#define HAVE_BLAS 1
#define HAVE_LAPACK 1
#define HAVE_LIBSTDC__ 1
#define HAVE_LIBM 1
#define HAVE_LIBPTHREAD 1
#define USING_GFORTRAN 1
#define NDEBUG 1
#define HAVE_MPI 1
#define HAVE_MPI_CXX 1
#define _MPI_CPP_BINDINGS 1
#define HAVE_LIBNETCDF 1
#define HAVE_NETCDF 1
#define HAVE_LIBNETCDFF 1
#define HAVE_PARMETIS 1
#define HAVE_PETSC_MODULES 1
#define HAVE_PETSC 1
#define HAVE_HYPRE 1
#if defined(HAVE_LIBUDUNITS) && defined(HAVE_LIBNETCDF) && defined(HAVE_LIBCGNS)
#define ENABLE_ERA40 1
#else
#error ERROR: ERA-40 support requested but support is missing for one of NetCDF, CGNS or Udunits
#endif
#define DOUBLEP 1
#define SIGNAL
#undef __FEMDEM_VERSION__
#define __FEMDEM_VERSION__ "Unversioned directory"
