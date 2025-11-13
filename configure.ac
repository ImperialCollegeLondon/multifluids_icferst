dnl -*- Autoconf -*-
dnl ICFERST/Fluidity modernized configure.ac (Ubuntu 24.04)
dnl Requires: autoconf>=2.69, automake>=1.16, libtool, autoconf-archive (for AX_*),
dnl           pkg-config, cmake
dnl   sudo apt install autoconf automake libtool autoconf-archive pkg-config cmake
dnl   sudo apt install petsc-dev libpetsc-real3.19-dev libhdf5-openmpi-dev libnetcdf-mpi-dev
dnl   sudo apt install libparmetis-dev libmetis-dev libexodusii-dev libudunits2-dev
dnl   sudo apt install libvtk9-dev python3-dev python3-numpy

AC_PREREQ([2.69])
AC_INIT([ICFERST],[0.0.0],[bugs@example.org])
AC_CONFIG_SRCDIR([main.cpp])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])

dnl Generate config header that your code includes as "include/config.h"
AC_CONFIG_HEADERS([include/config.h])

#AC_PROG_MAKE_SET
AC_PROG_INSTALL
AC_CHECK_PROG(MAKE, gmake, gmake, make, $PATH)

dnl ---- Language/toolchain basics
AC_PROG_CC
AC_PROG_CXX
AC_PROG_FC
AC_PROG_LN_S
AC_CANONICAL_HOST

if test "$FC" = "gfortran"; then
    AC_DEFINE([USING_GFORTRAN],[1],[Define to 1 if gfortran is the compiler])
fi

# Get the current date and Git commit
echo "character(len=*), parameter :: compile_date = '$(date +"%Y-%m-%d %H:%M:%S")'" > version.inc
echo "character(len=*), parameter :: git_commit = '$(git log -1 --pretty=format:"%H")'" >> version.inc

dnl ---- Let users override wrappers from env (MPIF90, MPIF77, MPICC, MPICXX)
AC_ARG_VAR([MPIF90], [MPI Fortran 90/2003 wrapper, e.g. mpifort/mpif90])
AC_ARG_VAR([MPIF77], [MPI Fortran 77 wrapper, e.g. mpif77])
AC_ARG_VAR([MPICC],  [MPI C wrapper, e.g. mpicc])
AC_ARG_VAR([MPICXX], [MPI C++ wrapper, e.g. mpicxx])

dnl Prefer wrappers if present; otherwise fall back to the base compilers
AC_PATH_PROGS([MPIF90], [mpifort mpif90 mpf90], [$FC])
AC_PATH_PROGS([MPIF77], [mpif77  mpf77],        [$F77])
AC_PATH_PROGS([MPICC],  [mpicc   mpcc],         [$CC])
AC_PATH_PROGS([MPICXX], [mpicxx  mpiCC mpCC],   [$CXX])

dnl Use the MPI wrappers for the actual build (like the old configure did)
F77="$MPIF77"
FC="$MPIF90"
CC="$MPICC"
CXX="$MPICXX"

AC_SUBST([F77])  AC_SUBST([FC])  AC_SUBST([CC])  AC_SUBST([CXX])
AC_SUBST([MPIF90]) AC_SUBST([MPIF77]) AC_SUBST([MPICC]) AC_SUBST([MPICXX])

dnl ---- Check we can link with the wrappers
AC_MSG_CHECKING([whether we can compile/link an MPI C program with $MPICC])
save_CC=$CC
CC=$MPICC
AC_LANG_PUSH([C])
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([[#include <mpi.h>]], [[int f=0; MPI_Initialized(&f); return 0;]])],
  [have_mpi_c=yes], [have_mpi_c=no])
AC_LANG_POP([C])
CC=$save_CC
AC_MSG_RESULT([$have_mpi_c])

AC_MSG_CHECKING([whether we can compile/link an MPI Fortran program with $MPIF90])
save_FC=$FC
FC=$MPIF90
AC_LANG_PUSH([Fortran])
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([], [[
      use mpi
      integer :: ierr, r
      call MPI_Init(ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD, r, ierr)
      call MPI_Finalize(ierr)
  ]])],
  [have_mpi_f=yes], [have_mpi_f=no])
AC_LANG_POP([Fortran])
FC=$save_FC
AC_MSG_RESULT([$have_mpi_f])

dnl Fallback path like your old configure: try -lmpi -lmpi_mpifh if wrappers fail
AS_IF([test "x$have_mpi_c" != xyes -o "x$have_mpi_f" != xyes], [
  AC_MSG_NOTICE([MPI wrappers failed; trying -lmpi fallback])
  save_LIBS=$LIBS
  LIBS="-lmpi -lmpi_mpifh $LIBS"
  AC_LANG_PUSH([C])
  save_CC=$CC; CC=${CC:-cc}
  AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([[#include <mpi.h>]], [[int f=0; MPI_Initialized(&f); return 0;]])],
    [have_mpi_c=yes], [have_mpi_c=no])
  CC=$save_CC
  AC_LANG_POP([C])
  LIBS=$save_LIBS
])

AS_IF([test "x$have_mpi_c" = xyes -a "x$have_mpi_f" = xyes], [
  AC_DEFINE([HAVE_MPI], [1], [Define if MPI is available])
], [
  AC_MSG_FAILURE([no MPI support found on system])
])

dnl ---- Optional: detect mpi_f08 (communicators as derived type)
AC_MSG_CHECKING([for Fortran mpi_f08 module])
save_FC=$FC
FC=$MPIF90
AC_LANG_PUSH([Fortran])
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM([], [[
      use mpi_f08
      type(MPI_Comm) :: c
      c = MPI_COMM_WORLD
  ]])],
  [have_mpi_f08=yes], [have_mpi_f08=no])
AC_LANG_POP([Fortran])
FC=$save_FC
AC_MSG_RESULT([$have_mpi_f08])
AS_IF([test "x$have_mpi_f08" = xyes], [
  AC_DEFINE([HAVE_MPI_F08], [1], [Define if mpi_f08 module is available])
  FCFLAGS="$FCFLAGS -DHAVE_MPI_F08"
])

dnl ---- Define Fortran/C name-mangling macros for your fmangle.h usage
dnl This ensures FC_FUNC / FC_FUNC_ land in config.h (map to F77_FUNC in your code)
AC_F77_WRAPPERS
AC_FC_FUNC([f77_dummy],[F77_DUMMY])


dnl -------- Libtool / ar --------
LT_INIT

dnl -------- Languages / tools --------
AC_PROG_CC
AC_PROG_CXX
AC_PROG_FC

dnl === gfortran >= 10: relax legacy interface arg checking (MPI wrappers etc.) ===
AC_CACHE_VAL([ac_cv_fc_version], [ac_cv_fc_version=`$FC -dumpfullversion -dumpversion 2>/dev/null`])
AC_MSG_NOTICE([Fortran compiler version: $ac_cv_fc_version])
FC_MAJOR=`echo "$ac_cv_fc_version" | sed 's/\..*$//'`

dnl Prefer feature test, but if not available just gate by version
AS_IF([test -n "$FC_MAJOR" -a "$FC_MAJOR" -ge 10], [
  dnl If autoconf-archive is present, verify the flag compiles; otherwise just add it.
  m4_ifdef([AX_CHECK_COMPILE_FLAG], [
    AX_CHECK_COMPILE_FLAG([-fallow-argument-mismatch],
      [FCFLAGS="$FCFLAGS -fallow-argument-mismatch"],
      [FCFLAGS="$FCFLAGS -fallow-argument-mismatch"])
  ], [
    FCFLAGS="$FCFLAGS -fallow-argument-mismatch"
  ])
])


dnl --- OpenMP detection (compile + link) ---
m4_ifdef([AX_OPENMP], [
  AX_OPENMP
  test -n "$OPENMP_CFLAGS"   && CFLAGS="$CFLAGS $OPENMP_CFLAGS"
  test -n "$OPENMP_CXXFLAGS" && CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
  test -n "$OPENMP_FCFLAGS"  && FCFLAGS="$FCFLAGS $OPENMP_FCFLAGS"
  dnl IMPORTANT: add to link too
  test -n "$OPENMP_LIBS"     && LIBS="$LIBS $OPENMP_LIBS"
  test -z "$OPENMP_LIBS" -a -n "$OPENMP_FCFLAGS" && LDFLAGS="$LDFLAGS $OPENMP_FCFLAGS"
], [
  dnl Fallback for GCC-like compilers
  CFLAGS="$CFLAGS -fopenmp"
  CXXFLAGS="$CXXFLAGS -fopenmp"
  FCFLAGS="$FCFLAGS -fopenmp"
  LDFLAGS="$LDFLAGS -fopenmp"
])

dnl --- Make sure these get substituted into Makefiles
AC_SUBST([OPENMP_FCFLAGS])
AC_SUBST([OPENMP_CFLAGS])
AC_SUBST([OPENMP_CXXFLAGS])
AC_SUBST([OPENMP_LIBS])
dnl Some AX_OPENMP variants don’t set OPENMP_LIBS; ensure a link flag is exported.
test -n "$OPENMP_CFLAGS" && OPENMP_LDFLAGS="$OPENMP_CFLAGS"
AC_SUBST([OPENMP_LDFLAGS])
LDFLAGS="$LDFLAGS $OPENMP_LDFLAGS"

dnl -- (optional) add them to your global flags if you want
FCFLAGS="$FCFLAGS $OPENMP_FCFLAGS"
CFLAGS="$CFLAGS $OPENMP_CFLAGS"
CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
LIBS="$LIBS $OPENMP_LIBS"

dnl --- We link C++ objects into the Fortran main; ensure C++ stdlib is on LIBS
AC_LANG_PUSH([C++])
AC_CHECK_LIB([stdc++], [__cxa_begin_catch],
  [LIBS="$LIBS -lstdc++"])
AC_LANG_POP([C++])

dnl ---- Disable legacy MPI C++ bindings to avoid REAL macro clashes ----
CPPFLAGS="$CPPFLAGS -DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX"

dnl --- Fortran module flags (add this block here) -------------------------
dnl Where to SEARCH for modules (include dirs):
MODINC_FLAG="-I"
AC_SUBST([MODINC_FLAG])

dnl Where to WRITE .mod files (compiler-specific):
AC_CACHE_CHECK([Fortran module output flag], [ac_cv_fc_modout_flag], [
  ac_cv_fc_modout_flag=
  AC_LANG_PUSH([Fortran])
  mkdir -p conftestmods
  save_FCFLAGS="$FCFLAGS"

  dnl gfortran-style: -J<dir>
  FCFLAGS="$save_FCFLAGS -Jconftestmods"
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[])],
                    [ac_cv_fc_modout_flag="-J"], [])

  dnl Intel/NVHPC-style: -module <dir>
  if test -z "$ac_cv_fc_modout_flag"; then
    FCFLAGS="$save_FCFLAGS -module conftestmods"
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[])],
                      [ac_cv_fc_modout_flag="-module "], [])
  fi

  FCFLAGS="$save_FCFLAGS"
  rm -rf conftestmods
  AC_LANG_POP([Fortran])

  test -z "$ac_cv_fc_modout_flag" && ac_cv_fc_modout_flag="-J"
])
MOD_FLAG="$ac_cv_fc_modout_flag"
AC_SUBST([MOD_FLAG])
dnl ------------------------------------------------------------------------
AC_LANG([C])
AC_PROG_INSTALL
dnl ---- Find archiver and ranlib (portable, no Automake required) ----
AC_CHECK_TOOL([AR],     [ar])
AC_CHECK_TOOL([RANLIB], [ranlib], [:])

dnl --- Archiver setup (static libraries) ---
AC_PATH_PROG([AR],[ar],[ar])
test -n "$ARFLAGS" || ARFLAGS="rcs"  dnl r=replace, c=create, s=index
AC_SUBST([AR])
AC_SUBST([ARFLAGS])
AC_PROG_SED
AC_PROG_GREP
PKG_PROG_PKG_CONFIG
dnl --- Optional loop-unrolling flag (used in some subdir Makefiles) ---
AC_ARG_ENABLE([unroll-loops],
  [AS_HELP_STRING([--enable-unroll-loops],
                  [Enable explicit loop-unroll flags (default: auto: gfortran -> -funroll-loops, others empty)])],
  [], [enable_unroll_loops=auto])

UNROLL_LOOPS=""

if test "x$enable_unroll_loops" = "xyes" -o "x$enable_unroll_loops" = "xauto"; then
  fc_version=$($FC --version 2>/dev/null | head -n1)
  case "$fc_version" in
    *gfortran*|*GNU*) UNROLL_LOOPS="-funroll-loops" ;;
    *ifort*|*Intel*)  UNROLL_LOOPS="-unroll" ;;
    # add other compilers here if you really want, otherwise leave empty
  esac
fi

AC_SUBST([UNROLL_LOOPS])

dnl -------- Options --------
AC_ARG_ENABLE([debugging],
  [AS_HELP_STRING([--enable-debugging],[Build with debug flags])],[],[enable_debugging=no])

AC_ARG_ENABLE([openmp],
  [AS_HELP_STRING([--enable-openmp],[Enable OpenMP])],[],[enable_openmp=no])

AC_ARG_ENABLE([profiling],
  [AS_HELP_STRING([--enable-profiling],[Enable gprof profiling (-pg)])],[],[enable_profiling=no])

AC_ARG_ENABLE([dp],
  [AS_HELP_STRING([--enable-dp[=FLAG]],[Force default REAL=8; if FLAG provided (e.g. -fdefault-real-8) it will be used])],[],[enable_dp=yes])

AC_ARG_ENABLE([2d-adaptivity],
  [AS_HELP_STRING([--enable-2d-adaptivity],[Link against 2D adaptivity library (GPL)])],[],[enable_2d_adaptivity=no])

AC_ARG_ENABLE([mba3d],
  [AS_HELP_STRING([--enable-mba3d],[Enable libmba3d 3D adaptivity library])],[],[enable_mba3d=no])

AC_ARG_ENABLE([memory-stats],
  [AS_HELP_STRING([--enable-memory-stats],[Enable memory statistics])],[],[enable_memory_stats=no])

AC_ARG_ENABLE([hyperlight],
  [AS_HELP_STRING([--enable-hyperlight],[Enable Hyperlight solar irradiance model])],[],[enable_hyperlight=no])

AC_ARG_ENABLE([vtk],
  [AS_HELP_STRING([--enable-vtk],[Enable VTK (default yes)])],[],[enable_vtk=yes])

AC_ARG_ENABLE([sam],
  [AS_HELP_STRING([--enable-sam],[Use SAM rather than Zoltan for adaptive load rebalancing])],[],[enable_sam=no])

dnl ---- Python embedding (robust, no stray semicolons) ----
AC_ARG_ENABLE([python],
  [AS_HELP_STRING([--enable-python],[embed Python in ICFERST (default: yes)])],
  [], [enable_python=yes])

AC_CHECK_PROGS([PYTHON],[python3 python], [python3])

AS_IF([test "x$enable_python" != "xno"], [

  dnl Prefer autoconf-archive’s AX_PYTHON_DEVEL if available
  m4_ifdef([AX_PYTHON_DEVEL], [
    AX_PYTHON_DEVEL([>= '3.8'])
    HAVE_PYTHON=yes
    AC_DEFINE([HAVE_PYTHON],[1],[Define to 1 to embed Python])

    CPPFLAGS="$CPPFLAGS $PYTHON_CPPFLAGS"
    CFLAGS="$CFLAGS $PYTHON_CPPFLAGS"
    CXXFLAGS="$CXXFLAGS $PYTHON_CPPFLAGS"
    LIBS="$LIBS $PYTHON_LIBS $PYTHON_EXTRA_LDFLAGS $PYTHON_EXTRA_LIBS"

  ], [
    dnl Fallback: use pkg-config (Ubuntu 24.04 provides python3-embed)
    m4_ifdef([PKG_CHECK_MODULES], [
      AC_MSG_NOTICE([AX_PYTHON_DEVEL not found; trying pkg-config python3-embed])
      PKG_CHECK_MODULES([PY3EMB], [python3-embed], [
        HAVE_PYTHON=yes
        AC_DEFINE([HAVE_PYTHON],[1],[Define to 1 to embed Python])
        CPPFLAGS="$CPPFLAGS $PY3EMB_CFLAGS"
        CFLAGS="$CFLAGS $PY3EMB_CFLAGS"
        CXXFLAGS="$CXXFLAGS $PY3EMB_CFLAGS"
        LIBS="$LIBS $PY3EMB_LIBS"
      ], [
        dnl Try versioned embed name as a fallback (e.g. python-3.12-embed)
        PKG_CHECK_MODULES([PY3EMB], [python-3.12-embed], [
          HAVE_PYTHON=yes
          AC_DEFINE([HAVE_PYTHON],[1],[Define to 1 to embed Python])
          CPPFLAGS="$CPPFLAGS $PY3_EMB_CFLAGS"
          CFLAGS="$CFLAGS $PY3_EMB_CFLAGS"
          CXXFLAGS="$CXXFLAGS $PY3_EMB_CFLAGS"
          LIBS="$LIBS $PY3_EMB_LIBS"
        ], [
          AC_MSG_WARN([No python embed development files found; disabling Python embedding])
          HAVE_PYTHON=no
        ])
      ])
      dnl Define the module path in config.h (as a C string) instead of injecting via -D into CFLAGS
      AC_DEFINE_UNQUOTED([FLUIDITY_PYTHON_MODULE_PATH], ["$(pwd)/python"], [Path to Fluidity Python module directory])
      AC_SUBST([FLUIDITY_PYTHON_MODULE_PATH])
    ], [
      AC_MSG_WARN([No AX_PYTHON_DEVEL or pkg-config macros available; disabling Python embedding])
      HAVE_PYTHON=no
    ])
  ])

  dnl (Optional) NumPy headers, non-fatal
  AS_IF([test "x$HAVE_PYTHON" = "xyes"], [
    NUMPY_INCL=`"$PYTHON" - <<'PY'
try:
  import numpy; print(numpy.get_include())
except Exception:
  pass
PY`
    AS_IF([test -n "$NUMPY_INCL"], [
      AC_MSG_NOTICE([NumPy headers at $NUMPY_INCL])
      CPPFLAGS="$CPPFLAGS -DHAVE_NUMPY -I$NUMPY_INCL"
      CFLAGS="$CFLAGS -DHAVE_NUMPY -I$NUMPY_INCL"
      CXXFLAGS="$CXXFLAGS -DHAVE_NUMPY -I$NUMPY_INCL"
      FFLAGS="$FFLAGS -DHAVE_NUMPY -I$NUMPY_INCL"
      FCFLAGS="$FCFLAGS -DHAVE_NUMPY -I$NUMPY_INCL"
    ], [
      AC_MSG_NOTICE([NumPy not found (optional)])
    ])
  ])

  dnl --- existing code ---

  dnl Define the module path in config.h (as a C string) instead of injecting via -D into CFLAGS
  AC_DEFINE_UNQUOTED([FLUIDITY_PYTHON_MODULE_PATH], ["$(pwd)/python"], [Path to Fluidity Python module directory])
  AC_SUBST([FLUIDITY_PYTHON_MODULE_PATH])

  dnl --- existing code ---
], [
  HAVE_PYTHON=no
])

AC_SUBST([HAVE_PYTHON])
dnl ---- end Python embedding ----

dnl Paths for optional components
#AC_ARG_WITH([vtk],
#  [AS_HELP_STRING([--with-vtk=DIR],[Path to VTKConfig.cmake (e.g. /usr/lib/x86_64-linux-gnu/cmake/vtk-9.1)])])

AC_ARG_WITH([zoltan],
  [AS_HELP_STRING([--with-zoltan=DIR],[Trilinos/Zoltan install prefix containing include/ and lib/])])

AC_ARG_WITH([exodusii],
  [AS_HELP_STRING([--with-exodusii@<:@=ROOT@:>@],[exodusII root (with include/ and lib/)])])

AC_ARG_WITH([judy],
  [AS_HELP_STRING([--with-judy=<lib>],[Use external Judy, e.g. -ljudy or /path/libjudy.a])])

AC_ARG_WITH([libspud-root],
  [AS_HELP_STRING([--with-libspud-root=<dir>],[Use external libspud root (with include/ and lib/)])])

AC_ARG_WITH([spatialindex-root],
  [AS_HELP_STRING([--with-spatialindex-root=<dir>],[Use external libspatialindex root])])

AC_ARG_WITH([cgal],
  [AS_HELP_STRING([--with-cgal],[Enable CGAL via pkg-config])])

AC_ARG_WITH([phreeqc],
  [AS_HELP_STRING([--with-phreeqc=DIR],[PHREEQC-CRM root, with lib/ and include/])])

AC_ARG_WITH([xgboost],
  [AS_HELP_STRING([--with-xgboost=DIR],[XGBoost root, with lib/ and include/])])

AC_ARG_WITH([psmile],
  [AS_HELP_STRING([--with-psmile=DIR],[OASIS/psmile root, with lib/ and include/])])

AC_ARG_WITH([libnuma],
  [AS_HELP_STRING([--with-libnuma='<LDFLAGS>'],[e.g. '--with-libnuma=-L/usr/lib -lnuma'])])

AC_ARG_WITH([femdem],
  [AS_HELP_STRING([--with-femdem=DIR],[FEMDEM root, with lib/])])

dnl -------- Host OS note --------
AC_CANONICAL_HOST
case $host_os in
  darwin*) AC_DEFINE([DARWIN],[1],[Building on Darwin]);;
  linux*)  AC_DEFINE([LINUX],[1],[Building on Linux]);;
esac

dnl -------- Feature-tested flags (needs autoconf-archive) --------
AX_CHECK_COMPILE_FLAG([-ffree-line-length-none],[FCFLAGS="$FCFLAGS -ffree-line-length-none"],[])
AX_CHECK_COMPILE_FLAG([-ffixed-line-length-none],[FCFLAGS="$FCFLAGS -ffixed-line-length-none"],[])

AS_IF([test "x$enable_debugging" = "xyes"], [
  AX_CHECK_COMPILE_FLAG([-Og],[CFLAGS="$CFLAGS -Og"; CXXFLAGS="$CXXFLAGS -Og"],[CFLAGS="$CFLAGS -O0"; CXXFLAGS="$CXXFLAGS -O0"])
  FCFLAGS="$FCFLAGS -g"
  CFLAGS="$CFLAGS -g"
  CXXFLAGS="$CXXFLAGS -g"
  AC_DEFINE([DDEBUG],[1],[Enable debug code])
], [
  AX_CHECK_COMPILE_FLAG([-O3],[CFLAGS="$CFLAGS -O3"; CXXFLAGS="$CXXFLAGS -O3"],[])
  AX_CHECK_COMPILE_FLAG([-O3],[FCFLAGS="$FCFLAGS -O3"],[])
  AC_DEFINE([NDEBUG],[1],[Disable asserts])
])

AS_IF([test "x$enable_profiling" = "xyes"], [
  FCFLAGS="-pg $FCFLAGS"
  CFLAGS="-pg $CFLAGS"
  CXXFLAGS="-pg $CXXFLAGS"
  LDFLAGS="-pg $LDFLAGS"
  AC_SUBST([PROFILING_FLAG],["-pg"])
], [AC_SUBST([PROFILING_FLAG],[""])])

AS_IF([test "x$enable_openmp" = "xyes"], [
  AX_OPENMP([],[AC_MSG_WARN([OpenMP not detected]); :])
  test -n "$OPENMP_CFLAGS" && CFLAGS="$CFLAGS $OPENMP_CFLAGS"
  test -n "$OPENMP_CXXFLAGS" && CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
  test -n "$OPENMP_FCFLAGS" && FCFLAGS="$FCFLAGS $OPENMP_FCFLAGS"
  AX_CHECK_COMPILE_FLAG([-fopenmp],[FCFLAGS="$FCFLAGS -fopenmp"],[])
])

dnl -------- Optional PIC toggle akin to old --enable-shared --------
AC_ARG_ENABLE([shared],
  [AS_HELP_STRING([--enable-shared],[Compile objects with -fPIC for shared builds])],[],[enable_shared=yes])
PIC_FLAG=""
AS_IF([test "x$enable_shared" = "xyes"],[
  AX_CHECK_COMPILE_FLAG([-fPIC],[PIC_FLAG="-fPIC"; CFLAGS="$CFLAGS -fPIC"; CXXFLAGS="$CXXFLAGS -fPIC"; FCFLAGS="$FCFLAGS -fPIC"],[])
])
AC_SUBST([PIC_FLAG])

dnl -------- Save pre-MPI compilers for subprojects that reconfigure --------
saved_CC="$CC"
saved_CXX="$CXX"
saved_FC="$FC"
F77="$FC"  dnl compat var for old sub-configures
AC_SUBST([F77])

dnl -------- MPI wrappers --------
AC_ARG_VAR([MPICC],[MPI C compiler])
AC_ARG_VAR([MPICXX],[MPI C++ compiler])
AC_ARG_VAR([MPIFC],[MPI Fortran compiler])

AC_PATH_PROGS([MPICC],[mpicc],[$CC])
AC_PATH_PROGS([MPICXX],[mpicxx],[$CXX])
AC_PATH_PROGS([MPIFC],[mpif90 mpifort],[$FC])

dnl --- Ensure we link helper libs for Fortran (mpif.h and/or mpi_f08) ---
AC_MSG_CHECKING([for Open MPI Fortran helper libs])
AC_CHECK_LIB([mpi_usempif08], [mpi_register_datarep_f08],
             [OMPI_HAVE_USEMPIF08=yes], [OMPI_HAVE_USEMPIF08=no])
AC_CHECK_LIB([mpi_mpifh],     [mpi_fortran_argv_null],
             [OMPI_HAVE_MPIFH=yes],     [OMPI_HAVE_MPIFH=no])
AC_MSG_RESULT([usempif08=$OMPI_HAVE_USEMPIF08 mpifh=$OMPI_HAVE_MPIFH])

MPI_HELPER_LIBS=""

dnl Force-keep DSOs even with --as-needed by referencing a symbol
AS_IF([test "x$OMPI_HAVE_USEMPIF08" = "xyes"], [
  MPI_HELPER_LIBS="$MPI_HELPER_LIBS -Wl,--no-as-needed -Wl,-u,mpi_register_datarep_f08 -lmpi_usempif08 -Wl,--as-needed"
])
AS_IF([test "x$OMPI_HAVE_MPIFH" = "xyes"], [
  MPI_HELPER_LIBS="$MPI_HELPER_LIBS -Wl,--no-as-needed -lmpi_mpifh -Wl,--as-needed"
])

AC_SUBST([MPI_HELPER_LIBS])

AS_IF([test -z "$MPICC" -o -z "$MPIFC"], [AC_MSG_ERROR([MPI compilers not found (need mpicc and mpif90/mpifort)])])

CC="$MPICC"
CXX="$MPICXX"
FC="$MPIFC"

AC_DEFINE([HAVE_MPI],[1],[Define if MPI is available])
AC_DEFINE([HAVE_MPI_CXX],[1],[Define if legacy MPI extensions are available])
AC_DEFINE([_MPI_CPP_BINDINGS],[1],[Define if legacy MPI extensions are available])
FLLINKER="$FC"
LINKER="$CXX"
AC_SUBST([FLLINKER])
AC_SUBST([LINKER])

AC_MSG_CHECKING([for MPI compile and link flags using mpicxx])
AS_IF([command -v mpicxx >/dev/null], [
    MPI_CXXFLAGS=`mpicxx --showme:compile`
    MPI_LIBS=`mpicxx --showme:link`
    AC_MSG_NOTICE([DEBUG: MPI flags: $MPI_CXXFLAGS, $MPI_LIBS])
], [
    MPI_CXXFLAGS=""
    MPI_LIBS=""
])
AC_SUBST([MPI_CXXFLAGS])
AC_SUBST([MPI_LIBS])

CPPFLAGS="$CPPFLAGS $MPI_CXXFLAGS"
LIBS="$LIBS $MPI_LIBS"

dnl -------- BLAS/LAPACK --------
have_blaslapack=no
PKG_CHECK_MODULES([BLAS],[blas],[
  CPPFLAGS="$CPPFLAGS $BLAS_CFLAGS"
  LIBS="$LIBS $BLAS_LIBS"
  have_blaslapack=yes
],[ : ])
PKG_CHECK_MODULES([LAPACK],[lapack],[
  CPPFLAGS="$CPPFLAGS $LAPACK_CFLAGS"
  LIBS="$LIBS $LAPACK_LIBS"
  have_blaslapack=yes
],[ : ])
AS_IF([test "x$have_blaslapack" = "xno"], [
  AC_CHECK_LIB([lapack],[dgesv],[],[AC_MSG_ERROR([LAPACK not found])])
  AC_CHECK_LIB([blas],[dgemm],[],[AC_MSG_ERROR([BLAS not found])])
])

dnl -------- HDF5 (prefer parallel) --------
PKG_CHECK_MODULES([HDF5],[hdf5-openmpi],[
  CPPFLAGS="$CPPFLAGS $HDF5_CFLAGS"
  LIBS="$LIBS $HDF5_LIBS"
],[
  PKG_CHECK_MODULES([HDF5],[hdf5],[
    CPPFLAGS="$CPPFLAGS $HDF5_CFLAGS"
    LIBS="$LIBS $HDF5_LIBS"
  ],[AC_MSG_WARN([HDF5 not found])])
])

dnl -------- NetCDF (C then Fortran) --------
PKG_CHECK_MODULES([NETCDF_C],[netcdf-mpi],[
  CPPFLAGS="$CPPFLAGS $NETCDF_C_CFLAGS"
  LIBS="$LIBS $NETCDF_C_LIBS"
],[
  PKG_CHECK_MODULES([NETCDF_C],[netcdf],[
    CPPFLAGS="$CPPFLAGS $NETCDF_C_CFLAGS"
    LIBS="$LIBS $NETCDF_C_LIBS"
  ],[AC_MSG_WARN([NetCDF C not found])])
])
PKG_CHECK_MODULES([NETCDFF],[netcdf-fortran],[
  CPPFLAGS="$CPPFLAGS $NETCDFF_CFLAGS"
  LIBS="$LIBS $NETCDFF_LIBS"
],[AC_MSG_WARN([NetCDF Fortran (netcdf-fortran) not found])])

dnl -------- UDUNITS2 (optional) --------
PKG_CHECK_MODULES([UDUNITS2],[udunits2],[
  CPPFLAGS="$CPPFLAGS $UDUNITS2_CFLAGS"
  LIBS="$LIBS $UDUNITS2_LIBS"
  AC_DEFINE([HAVE_LIBUDUNITS],[1],[Have UDUNITS2])
],[ : ])
AC_CHECK_LIB([udunits2],[utInit],[
  LIBS="$LIBS -ludunits2"
  AC_DEFINE([HAVE_LIBUDUNITS],[1],[Have UDUNITS2])
],[ : ])

PKG_PROG_PKG_CONFIG

dnl -------- PETSc --------
m4_ifdef([PKG_CHECK_MODULES], [
  AC_MSG_CHECKING([for PETSc via pkg-config])
  PKG_CHECK_MODULES([PETSC],[petsc >= 3.12],[
    CPPFLAGS="$CPPFLAGS $PETSC_CFLAGS"
    CFLAGS="$CFLAGS $PETSC_CFLAGS"
    CXXFLAGS="$CXXFLAGS $PETSC_CFLAGS"
    FCFLAGS="$FCFLAGS $PETSC_CFLAGS"
    LDFLAGS="$LDFLAGS $PETSC_LIBS"
    LIBS="$LIBS $PETSC_LIBS"
    AC_DEFINE([HAVE_PETSC],[1],[Have PETSc])
    HAVE_PETSC=yes
  ],[
    AC_MSG_WARN([pkg-config did not find a 'petsc' package])
    HAVE_PETSC=no
  ])
], [
  AC_MSG_WARN([PKG_CHECK_MODULES macro not available; will try PETSC_DIR fallback])
  HAVE_PETSC=no
])

dnl Fallback: respect PETSC_DIR if set
AS_IF([test "x$HAVE_PETSC" != "xyes" -a -n "$PETSC_DIR"], [
  AC_MSG_NOTICE([Using PETSC_DIR=$PETSC_DIR])
  CPPFLAGS="$CPPFLAGS -I$PETSC_DIR/include"
  CFLAGS="$CFLAGS -I$PETSC_DIR/include"
  CXXFLAGS="$CXXFLAGS -I$PETSC_DIR/include"
  FCFLAGS="$FCFLAGS -I$PETSC_DIR/include"
  LDFLAGS="$LDFLAGS -L$PETSC_DIR/lib"
  LIBS="$LIBS -L$PETSC_DIR/lib -lpetsc"
  AC_DEFINE([HAVE_PETSC],[1],[Have PETSc])
  HAVE_PETSC=yes
])

AS_IF([test "x$HAVE_PETSC" = "xyes"], [
  AC_MSG_RESULT([PETSc: yes])
  AC_SUBST([HAVE_PETSC])
  AC_SUBST([PETSC_CFLAGS])
  AC_SUBST([PETSC_LIBS])
], [
  AC_MSG_ERROR([PETSc not found; install petsc-dev or set PETSC_DIR])
])

dnl -------- HYPRE (optional) --------
have_hypre=no
HYPRE_LIBS=""
HYPRE_CFLAGS=""

dnl Prefer pkg-config if available
m4_ifdef([PKG_CHECK_MODULES], [
  PKG_CHECK_MODULES([HYPRE],[hypre],[
    have_hypre=yes
  ],[
    dnl No .pc or not found: fallback to manual header/lib probe
    have_hypre=no
  ])
],[
  have_hypre=no
])

AS_IF([test "x$have_hypre" = "xno"], [
  dnl Common include dirs on Debian/Ubuntu
  for d in /usr/include/hypre /usr/include/x86_64-linux-gnu/hypre ; do
    if test -f "$d/HYPRE.h" ; then
      HYPRE_CFLAGS="$HYPRE_CFLAGS -I$d"
      break
    fi
  done

  dnl Try both include-names; one of these will be true if we set HYPRE_CFLAGS
  save_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS $HYPRE_CFLAGS"
  AC_CHECK_HEADERS([HYPRE.h hypre/HYPRE.h],[
    AC_CHECK_LIB([HYPRE],[HYPRE_IJMatrixCreate],[
      HYPRE_LIBS="$HYPRE_LIBS -lHYPRE"
      have_hypre=yes
    ],[])
  ],[])
  CPPFLAGS="$save_CPPFLAGS"
])

AS_IF([test "x$have_hypre" = "xyes"],[
  CPPFLAGS="$CPPFLAGS $HYPRE_CFLAGS"
  LIBS="$LIBS $HYPRE_LIBS"
  AC_DEFINE([HAVE_HYPRE],[1],[Have Hypre preconditioners])
])

AC_SUBST([HYPRE_CFLAGS])
AC_SUBST([HYPRE_LIBS])


dnl -------- METIS / ParMETIS --------
AC_CHECK_HEADERS([metis.h])
AC_CHECK_HEADERS([parmetis.h])
AC_CHECK_LIB([metis],[METIS_PartGraphKway],[AC_DEFINE([HAVE_LIBMETIS],[1],[Have METIS])],[])
AC_CHECK_LIB([parmetis],[ParMETIS_V3_AdaptiveRepart],[
  AC_DEFINE([HAVE_LIBPARMETIS],[1],[Have ParMETIS])
  LIBS="$LIBS -lparmetis -lmetis"
],[])
AC_LANG_PUSH([C])
AC_MSG_CHECKING([if ParMETIS major version <= 3])
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM([[
    #include <parmetis.h>
    #ifndef PARMETIS_MAJOR_VERSION
    #error "no version macro"
    #endif
    #if PARMETIS_MAJOR_VERSION > 3
    #error "ParMETIS >= 4"
    #endif
  ]], [[]])],
  [ AC_MSG_RESULT([yes])
    AC_DEFINE([HAVE_PARMETIS_V3],[1],[Define if ParMETIS major version is 3 or lower])
    PARMETIS_FLAGS="-DPARMETIS_V3"
  ],
  [ AC_MSG_RESULT([no])
    PARMETIS_FLAGS=""
  ])
AC_LANG_POP([C])
CPPFLAGS="$CPPFLAGS $PARMETIS_FLAGS"

dnl -------- ARPACK (optional) --------
PKG_CHECK_MODULES([ARPACK],[arpack],[
  CPPFLAGS="$CPPFLAGS $ARPACK_CFLAGS"
  LIBS="$LIBS $ARPACK_LIBS"
],[ AC_CHECK_LIB([arpack],[dsaupd],[],[ : ]) ])

dnl -------- Zoltan vs SAM --------
HAVE_ZOLTAN=no
AS_IF([test "x$enable_sam" = "xyes"], [
  HAVE_ZOLTAN=no
],[
  AS_IF([test -n "$with_zoltan"], [
    CPPFLAGS="$CPPFLAGS -I$with_zoltan/include"
    FCFLAGS="$FCFLAGS -I$with_zoltan/include"
    LDFLAGS="$LDFLAGS -L$with_zoltan/lib"
  ])
  AC_CHECK_HEADERS([zoltan.h],[],[have_zoltan_hdr=no])
  AC_MSG_CHECKING([whether to link against zoltan or trilinos_zoltan])
  have_zoltan_lib=no
  AC_CHECK_LIB([trilinos_zoltan],[Zoltan_Create],[have_zoltan_lib=yes; ZOLTAN_LIBS="-ltrilinos_zoltan"],[
    AC_CHECK_LIB([zoltan],[Zoltan_Create],[have_zoltan_lib=yes; ZOLTAN_LIBS="-lzoltan"],[
      AC_MSG_WARN([Neither zoltan nor trilinos_zoltan library found])
    ])
  ])
  AC_SUBST([ZOLTAN_LIBS])
  LIBS="$LIBS $ZOLTAN_LIBS"
  AS_IF([test "x$have_zoltan_hdr" != "xno" -a "x$have_zoltan_lib" != "xno"], [
    AC_LANG_PUSH([Fortran])
    AC_MSG_CHECKING([for usable zoltan Fortran module])
    cat > conftest_zoltan.F90 <<'EOT'
program test_zoltan
  use zoltan
end program test_zoltan
EOT
    ${FC} $FCFLAGS -c conftest_zoltan.F90 >/dev/null 2>&1
    if test "$?" -eq 0 ; then
      AC_MSG_RESULT([yes])
      AC_DEFINE([HAVE_ZOLTAN],[1],[Have Zoltan])
      HAVE_ZOLTAN=yes
    else
      AC_MSG_RESULT([no])
      AC_MSG_WARN([Zoltan library found but Fortran module unusable; rebuild Trilinos with same compiler/MPI])
    fi
    AC_LANG_POP([Fortran])
  ],[:])
])
AC_SUBST([HAVE_ZOLTAN])

dnl -------- CGAL (optional) --------
AS_IF([test "x$with_cgal" != "xno"],[
  PKG_CHECK_MODULES([CGAL],[cgal],[
    CPPFLAGS="$CPPFLAGS $CGAL_CFLAGS"
    LIBS="$LIBS $CGAL_LIBS -lgmp -lmpfr"
    AC_DEFINE([HAVE_LIBCGAL],[1],[Have CGAL])
  ],[
    AC_CHECK_LIB([CGAL],[_ZNK4CGAL6Random8get_seedEv],[LIBS="$LIBS -lCGAL -lgmp -lmpfr"; AC_DEFINE([HAVE_LIBCGAL],[1])],
      [AC_MSG_WARN([CGAL not found])])
  ])
],[:])

dnl -------- Exodus-II (optional) --------
HAVE_EXODUSII=no
AS_IF([test -n "$with_exodusii" -a "x$with_exodusii" != "xyes"],[
  CPPFLAGS="$CPPFLAGS -I$with_exodusii/include"
  LDFLAGS="$LDFLAGS -L$with_exodusii/lib"
],[:])
PKG_CHECK_MODULES([EXODUSII],[exodusii],[
  CPPFLAGS="$CPPFLAGS $EXODUSII_CFLAGS"
  LIBS="$LIBS $EXODUSII_LIBS"
  HAVE_EXODUSII=yes
  AC_DEFINE([HAVE_EXODUSII],[1],[Have Exodus-II])
],[ dnl fallback
  AC_CHECK_LIB([exodus],[ex_close],[LIBS="$LIBS -lexodus"; HAVE_EXODUSII=yes; AC_DEFINE([HAVE_EXODUSII],[1])],[ : ])
])
AC_SUBST([HAVE_EXODUSII])

dnl -------- PHREEQC-CRM (optional) --------
PHREEQC=""
AS_IF([test -n "$with_phreeqc"],[
  CPPFLAGS="$CPPFLAGS -I$with_phreeqc/include"
  LDFLAGS="$LDFLAGS -L$with_phreeqc/lib"
  LIB_PHREEQC="$with_phreeqc/lib/libphreeqcrm.so"
  PHREEQC="yes"
  AC_DEFINE([USING_PHREEQC],[1],[Use PHREEQC])
])
AC_SUBST([PHREEQC])
AC_SUBST([LIB_PHREEQC])

dnl -------- XGBOOST (optional) --------
XGBOOST=""
AS_IF([test -n "$with_xgboost"],[
  CPPFLAGS="$CPPFLAGS -I$with_xgboost/include"
  LDFLAGS="$LDFLAGS -L$with_xgboost/lib"
  LIB_XGBOOST="$with_xgboost/lib/libxgboost.so"
  XGBOOST="yes"
  AC_DEFINE([USING_XGBOOST],[1],[Use XGBoost])
])
AC_SUBST([XGBOOST])
AC_SUBST([LIB_XGBOOST])

dnl -------- OASIS/psmile (optional) --------
ENABLE_PSMILE=""
AS_IF([test -n "$with_psmile"],[
  CPPFLAGS="$CPPFLAGS -I$with_psmile/include"
  LDFLAGS="$LDFLAGS -L$with_psmile/lib"
  ENABLE_PSMILE="yes"
  AC_DEFINE([USING_PSMILE],[1],[Use OASIS psmile])
])
AC_SUBST([ENABLE_PSMILE])

dnl -------- FEMDEM (optional) --------
FEMDEM=""
AS_IF([test -n "$with_femdem"],[
  FLFEMDEM_PATH="$with_femdem"
  LIB_FEMDEM="$with_femdem/lib/libdfemdem.a"
  FEMDEM="yes"
  AC_DEFINE([USING_FEMDEM],[1],[Use FEMDEM])
])
AC_SUBST([FEMDEM])
AC_SUBST([FLFEMDEM_PATH])
AC_SUBST([LIB_FEMDEM])

dnl -------- Default REAL=8 (optional) --------
AS_IF([test "x$enable_dp" != "xno"], [
  if test "x$enable_dp" = "xyes"; then
    AX_CHECK_COMPILE_FLAG([-fdefault-real-8],[FCFLAGS="$FCFLAGS -fdefault-real-8 -fdefault-double-8"; AC_DEFINE([DOUBLEP],[1],[Default REAL=8])],[])
  else
    FCFLAGS="$FCFLAGS $enable_dp"
    AC_DEFINE([DOUBLEP],[1],[Default REAL=8])
  fi
])

dnl -------- Fortran stream I/O probe (harmless) --------
AC_LANG_PUSH([Fortran])
AC_MSG_CHECKING([Fortran stream I/O support])
AC_LINK_IFELSE([AC_LANG_SOURCE([
      program teststream
      integer :: u
      u=9
      open(unit=u, access="stream")
      end program teststream
])],[AC_MSG_RESULT([yes]); AC_DEFINE([STREAM_IO],[1],[Fortran stream I/O])],
  [AC_MSG_RESULT([no])])
AC_LANG_POP([Fortran])

dnl -------- Fortran signal probe (harmless) --------
AC_LANG_PUSH([Fortran])

AC_MSG_CHECKING([Fortran signal() support])

dnl Try the 3-argument variant first (has the extra flag)
AC_LINK_IFELSE(
  [AC_LANG_SOURCE([[
      program testsignal3
      implicit none
      integer :: r
      interface
        function handle_sigint(signum)
          integer :: handle_sigint
          integer, intent(in) :: signum
        end function handle_sigint
      end interface
      r = signal(2, handle_sigint, -1)
      end program testsignal3

      function handle_sigint(signum)
        implicit none
        integer :: handle_sigint
        integer, intent(in) :: signum
        handle_sigint = 0
      end function handle_sigint
  ]])],
  [ dnl success: 3-arg variant exists
    AC_MSG_RESULT([yes (3-arg)])
    AC_DEFINE([SIGNAL],[1],[Define to 1 if Fortran provides signal()])
    AC_DEFINE([SIGNAL_HAVE_FLAG],[1],[Define to 1 if signal() accepts a third flag argument])
  ],
  [ dnl otherwise try the 2-argument variant
    AC_LINK_IFELSE(
      [AC_LANG_SOURCE([[
          program testsignal2
          implicit none
          integer :: r
          interface
            function handle_sigint(signum)
              integer :: handle_sigint
              integer, intent(in) :: signum
            end function handle_sigint
          end interface
          r = signal(2, handle_sigint)
          end program testsignal2

          function handle_sigint(signum)
            implicit none
            integer :: handle_sigint
            integer, intent(in) :: signum
            handle_sigint = 0
          end function handle_sigint
      ]])],
      [ AC_MSG_RESULT([yes (2-arg)])
        AC_DEFINE([SIGNAL],[1],[Define to 1 if Fortran provides signal()])
        dnl leave SIGNAL_HAVE_FLAG undefined
      ],
      [ AC_MSG_RESULT([no]) ]
    )
  ]
)

AC_LANG_POP([Fortran])

dnl -------- Optional NUMA --------
AS_IF([test -n "$with_libnuma"],[
  AC_MSG_NOTICE([Linking against $with_libnuma])
  LIBS="$with_libnuma $LIBS"
  CPPFLAGS="-DHAVE_LIBNUMA $CPPFLAGS"
])

dnl -------- Optional tcmalloc --------
AC_CHECK_LIB([tcmalloc_minimal],[malloc],[LIBS="$LIBS -ltcmalloc_minimal"],[ : ])

dnl =========================
dnl  VTK via CMake (default on)
dnl =========================
AC_ARG_ENABLE([vtk],
  [AS_HELP_STRING([--disable-vtk],[Build without VTK (default: enabled)])],
  [enable_vtk="$enableval"], [enable_vtk=yes])

AC_ARG_WITH([vtk],
  [AS_HELP_STRING([--with-vtk=DIR],
                  [Path to VTKConfig.cmake (e.g. /usr/lib/x86_64-linux-gnu/cmake/vtk-9.1)])])

VTK_CFLAGS=""
VTK_LIBS=""

AS_IF([test "x$enable_vtk" != "xno"], [
  AC_PATH_PROG([CMAKE],[cmake],[])
  AS_IF([test -z "$CMAKE"], [AC_MSG_ERROR([cmake not found; install cmake])])

  dnl If no --with-vtk given, try common places for VTKConfig.cmake
  AS_IF([test -z "$with_vtk"], [
    for d in /usr/lib/x86_64-linux-gnu/cmake/vtk-* /usr/lib/cmake/vtk-* ; do
      test -f "$d/VTKConfig.cmake" -o -f "$d/vtk-config.cmake" && with_vtk="$d" && break
    done
  ])

  dnl Require a VTK config dir
  AS_IF([test ! -f "$with_vtk/VTKConfig.cmake" -a ! -f "$with_vtk/vtk-config.cmake"], [
    AC_MSG_WARN([VTKConfig.cmake not found at $with_vtk; will try header-only fallback])
  ])

  dnl Try CMake’s --find-package first
  AC_MSG_CHECKING([for VTK compile flags (via CMake)])
  VTK_CFLAGS=`$CMAKE -DNAME=VTK -DLANGUAGE=CXX -DMODE=COMPILE -DQUIET=1 \
                      -DCMAKE_PREFIX_PATH="$with_vtk" --find-package 2>/dev/null`; rc=$?
  AS_IF([test $rc -ne 0 || echo "$VTK_CFLAGS" | grep -qi 'not found'], [VTK_CFLAGS=""])
  AC_MSG_RESULT([${VTK_CFLAGS:-"(none from CMake)"}])

  AC_MSG_CHECKING([for VTK link flags (via CMake)])
  VTK_LIBS=`$CMAKE -DNAME=VTK -DLANGUAGE=CXX -DMODE=LINK -DQUIET=1 \
                   -DCMAKE_PREFIX_PATH="$with_vtk" --find-package 2>/dev/null`; rc=$?
  AS_IF([test $rc -ne 0 || echo "$VTK_LIBS" | grep -qi 'not found'], [VTK_LIBS=""])
  AC_MSG_RESULT([${VTK_LIBS:-"(none from CMake)"}])

  dnl -------- fallback include path for Ubuntu 24.04 VTK 9.1 --------
  AS_IF([test -z "$VTK_CFLAGS"], [
    for inc in /usr/include/vtk-9.1 /usr/local/include/vtk-9.1 ; do
      if test -d "$inc" ; then
        VTK_CFLAGS="-I$inc"
        break
      fi
    done
  ])

  dnl If still nothing, stop now (don’t leak text into flags)
  AS_IF([test -z "$VTK_CFLAGS"], [
    AC_MSG_ERROR([Could not determine VTK headers. Install libvtk9-dev or pass --with-vtk=/path/to/vtk-9.1 CMake dir])
  ])
  AS_IF([test -z "$VTK_LIBS"], [
    VTK_LIBS_ARRAY=""
    for lib in vtkCommonCore vtkCommonDataModel vtkCommonExecutionModel vtkFiltersCore vtkFiltersGeneral \
                vtkFiltersVerdict vtkIOCore vtkIOXML vtkIOGeometry vtkIOLegacy vtkIOParallelXML \
                vtkParallelMPI vtkIOXMLParser vtkCommonSystem vtkCommonMisc vtkCommonMath vtksys; do
      VTK_LIBS_ARRAY="$VTK_LIBS_ARRAY -l$lib-9.1"
    done
    VTK_LIBS="$VTK_LIBS_ARRAY"
  ])

  AC_SUBST([VTK_CFLAGS])
  AC_SUBST([VTK_LIBS])

  dnl NOTE: add include flags to all three; link flags only to LIBS
  CPPFLAGS="$CPPFLAGS $VTK_CFLAGS"
  CXXFLAGS="$CXXFLAGS $VTK_CFLAGS"
  FCFLAGS="$FCFLAGS $VTK_CFLAGS"
  LIBS="$LIBS $VTK_LIBS"

  AC_DEFINE([HAVE_VTK],[1],[Define to 1 if VTK is available])
])

AC_SUBST([VTK_CFLAGS])
AC_SUBST([VTK_LIBS])

dnl -------- Memory stats switches --------
AS_IF([test "x$enable_debugging" = "xyes"], [
  AS_IF([test "x$enable_memory_stats" != "xno"], [
    echo "#define HAVE_MEMORY_STATS 1" >> include/confdefs.h
    MEMORY_STATS="yes"
  ])
],[
  AS_IF([test "x$enable_memory_stats" = "xyes"], [
    echo "#define HAVE_MEMORY_STATS 1" >> include/confdefs.h
    MEMORY_STATS="yes"
  ])
])
AC_SUBST([MEMORY_STATS])

dnl -------- Hyperlight module --------
AS_IF([test "x$enable_hyperlight" = "xyes"],[
  AC_CHECK_FILES([hyperlight/Hyperlight.h hyperlight/Hyperlight.cpp hyperlight/Sky.h hyperlight/Sky.cpp hyperlight/Iop.h hyperlight/Iop.cpp],[
    echo "#define HAVE_HYPERLIGHT 1" >> include/confdefs.h
    HYPERLIGHT="yes"
    AC_SUBST([HYPERLIGHT])
  ],[
    AC_MSG_ERROR([Hyperlight enabled but sources missing.])
  ])
])

dnl -------- Subproject: libadaptivity --------
AS_IF([test -d libadaptivity],[
  AC_MSG_NOTICE([Configuring libadaptivity])
  dpflag=""
  test "x$enable_dp" = "xno" && dpflag="--enable-dp=no"
  ( cd libadaptivity && \
    FC="${saved_FC}" CC="${saved_CC}" CXX="${saved_CXX}" \
    FFLAGS="${PIC_FLAG} ${PROFILING_FLAG}" FCFLAGS="${PIC_FLAG} ${PROFILING_FLAG}" \
    CFLAGS="${PIC_FLAG} ${PROFILING_FLAG}"  CXXFLAGS="${PIC_FLAG} ${PROFILING_FLAG}" \
    CPPFLAGS="" ./configure ${dpflag} ) || AC_MSG_ERROR([libadaptivity configure failed])
  BUILD_LIBADAPT='@echo "    MAKE libadaptivity"; $(MAKE) -s -C libadaptivity && cp libadaptivity/lib/libadaptivity.a lib/'
  CLEAN_LIBADAPT='@echo "    CLEAN libadaptivity"; $(MAKE) -s -C libadaptivity clean'
  AC_SUBST([BUILD_LIBADAPT])
  AC_SUBST([CLEAN_LIBADAPT])
  AC_DEFINE([HAVE_ADAPTIVITY],[1],[Have libadaptivity])
],[
  AC_MSG_ERROR([libadaptivity directory missing])
])

dnl -------- Subproject: libjudy (internal or external) --------
FLIBJUDY=""
AS_IF([test -n "$with_judy"],[
  AC_MSG_NOTICE([Linking external Judy: $with_judy])
  case $with_judy in
    -l* | */* | *.a | *.so | *.so.* | *.o) acx_judy_lib="$with_judy" ;;
    *) acx_judy_lib="-l$with_judy" ;;
  esac
  FLIBJUDY="$acx_judy_lib"
],[
  AS_IF([test -d libjudy],[
    AC_MSG_NOTICE([Configuring libjudy])
    ( cd libjudy && \
      FC="${saved_FC}" CC="${saved_CC}" FFLAGS="${PIC_FLAG} ${PROFILING_FLAG}" \
      ./configure --prefix="$(pwd)/.." ) || AC_MSG_ERROR([libjudy configure failed])
    AC_LANG_PUSH([C])
    AC_CHECK_SIZEOF([void *])
    AC_LANG_POP([C])
    if test "$ac_cv_sizeof_void_p" = 8; then
      JUDY_CFLAGS="-DJU_64BIT"
    else
      JUDY_CFLAGS="-UJU_64BIT"
    fi
    echo "$CFLAGS" | grep -q -- "-O[23]" && JUDY_CFLAGS="$JUDY_CFLAGS -fno-aggressive-loop-optimizations"
    AC_SUBST([JUDY_CFLAGS])
    LIBJUDY=-ljudy
    FLIBJUDY="./lib/libjudy.a"
    AC_DEFINE([HAVE_JUDY],[1],[Have Judy])
  ],[
    AC_MSG_ERROR([libjudy directory missing and no --with-judy provided])
  ])
])
AC_SUBST([FLIBJUDY])

dnl -------- MBA2D/MBA3D optional libs --------
AS_IF([test "x$enable_2d_adaptivity" = "xyes"],[
  AC_DEFINE([HAVE_MBA_2D],[1],[2d adaptivity enabled])
  LIBS="-lmba2d $LIBS"
  MBA2D="yes"
  AC_SUBST([MBA2D])
])
AS_IF([test "x$enable_mba3d" = "xyes"],[
  AC_DEFINE([HAVE_MBA_3D],[1],[3d adaptivity enabled])
  LIBS="$LIBS -lmba3d"
  MBA3D="yes"
  AC_SUBST([MBA3D])
])

dnl -------- Subproject: libspud (internal or external) --------
AS_IF([test -n "$with_libspud_root"],[
  SPUDLIB="${with_libspud_root}/lib/libspud.a"
  SPUDHEADER="${with_libspud_root}/include/spud.h"
  CPPFLAGS="$CPPFLAGS -I${with_libspud_root}/include"
  FCFLAGS="$FCFLAGS -I${with_libspud_root}/include"
],[
  AS_IF([test -d libspud],[
    AC_MSG_NOTICE([Configuring libspud])
    CONFIG_FLAG="--disable-shared --enable-static"
    test "x$enable_profiling" = "xyes" && CONFIG_FLAG="${CONFIG_FLAG} --enable-profiling"
    SPUDLIB="${PWD}/lib/libspud.a"
    SPUDHEADER="${PWD}/include/spud.h"
    # ensure we’re not reusing an MPI-configured libspud
    ( cd libspud && \
      CC="${saved_CC}" CXX="${saved_CXX}" FC="${saved_FC}" \
      LIBS="" LDFLAGS="" \
      FFLAGS="${PIC_FLAG} ${PROFILING_FLAG}" FCFLAGS="${PIC_FLAG} ${PROFILING_FLAG}" \
      CFLAGS="${PIC_FLAG} ${PROFILING_FLAG}"  CXXFLAGS="${PIC_FLAG} ${PROFILING_FLAG}" \
      CPPFLAGS="" ./configure --prefix="$(pwd)/.." ${CONFIG_FLAG} ) \
      || AC_MSG_ERROR([libspud configure failed])
 
    dnl *** FORCE plain compilers on the sub-make as well ***
    BUILD_LIBSPUD='@echo "    MAKE libspud"; \
      $(MAKE) -s -C libspud \
        CC="$(saved_CC)" CXX="$(saved_CXX)" FC="$(saved_FC)" \
        CFLAGS="$(PIC_FLAG) $(PROFILING_FLAG)" \
        CXXFLAGS="$(PIC_FLAG) $(PROFILING_FLAG)" \
        FFLAGS="$(PIC_FLAG) $(PROFILING_FLAG)" \
        FCFLAGS="$(PIC_FLAG) $(PROFILING_FLAG)"; \
      $(MAKE) -s -C libspud install-libspud'

    CLEAN_LIBSPUD='@echo "    CLEAN libspud"; $(MAKE) -s -C libspud clean'
    AC_SUBST([BUILD_LIBSPUD])
    AC_SUBST([CLEAN_LIBSPUD])
    AC_DEFINE([HAVE_SPUD],[1],[Have libspud])
  ],[
    AC_MSG_ERROR([libspud directory missing and no --with-libspud-root provided])
  ])
])
AC_SUBST([SPUDLIB])
AC_SUBST([SPUDHEADER])

dnl -------- Subproject: spatialindex --------
AS_IF([test -n "$with_spatialindex_root"],[
  LIBSPATIALINDEX="${with_spatialindex_root}/lib/libspatialindex.a"
  CPPFLAGS="$CPPFLAGS -I${with_spatialindex_root}/include/spatialindex"
],[
  SPATIALINDEXDIR="spatialindex-1.8.0"
  AC_SUBST([SPATIALINDEXDIR])
  AS_IF([test -d "$SPATIALINDEXDIR"],[
    AC_MSG_NOTICE([Configuring spatialindex])
    ( cd "$SPATIALINDEXDIR" && \
      CC="${saved_CC}" CXX="${saved_CXX}" CFLAGS="" CXXFLAGS="" CPPFLAGS="-D_GNU_SOURCE" \
      ./configure --prefix="$(pwd)/../" --disable-shared --with-pic \
      --disable-maintainer-mode ) || AC_MSG_ERROR([spatialindex configure failed])
    LIBSPATIALINDEX="${PWD}/lib/libspatialindex.a"
  ],[
    AC_MSG_ERROR([spatialindex-1.8.0 directory missing; or pass --with-spatialindex-root])
  ])
])
AC_SUBST([LIBSPATIALINDEX])

dnl -------- Confdefs propagation like original --------
AS_IF([test ! -d include],[mkdir -p include])
grep "^#define DOUBLEP" confdefs.h > include/confdefs.h 2>/dev/null || true
grep "^#define F77_FUNC" confdefs.h >> include/confdefs.h 2>/dev/null || true
grep "^#define VTK" confdefs.h >> include/confdefs.h 2>/dev/null || true
grep "^#define MPI" confdefs.h >> include/confdefs.h 2>/dev/null || true
grep "^#define HAVE_HYPRE" confdefs.h >> include/confdefs.h 2>/dev/null || true
AS_IF([test ! -d libvtkfortran/include],[mkdir -p libvtkfortran/include])
grep "^#define DOUBLEP" confdefs.h > libvtkfortran/include/confdefs.h 2>/dev/null || true
grep "^#define F77_FUNC" confdefs.h >> libvtkfortran/include/confdefs.h 2>/dev/null || true
grep "^#define VTK" confdefs.h >> libvtkfortran/include/confdefs.h 2>/dev/null || true
grep "^#define MPI" confdefs.h >> libvtkfortran/include/confdefs.h 2>/dev/null || true
grep "^#define HAVE_HYPRE" confdefs.h >> libvtkfortran/include/confdefs.h 2>/dev/null || true

echo "#define USING_CONFIG_H 1" >> include/confdefs.h
echo "#include \"config.h\"" >> include/confdefs.h
echo "#define USING_CONFIG_H 1" >> libvtkfortran/include/confdefs.h
echo "#include \"config.h\"" >> libvtkfortran/include/confdefs.h

dnl -------- Summary --------
AC_MSG_NOTICE([==== Summary ====])
AC_MSG_NOTICE([CC=$CC])
AC_MSG_NOTICE([CXX=$CXX])
AC_MSG_NOTICE([FC=$FC])
AC_MSG_NOTICE([CFLAGS=$CFLAGS])
AC_MSG_NOTICE([CXXFLAGS=$CXXFLAGS])
AC_MSG_NOTICE([FCFLAGS=$FCFLAGS])
AC_MSG_NOTICE([LDFLAGS=$LDFLAGS])
AC_MSG_NOTICE([PETSc: yes])
dnl ----- Summary booleans (avoid $(test ... ) inside AC_MSG_*) -----
hdf5_have=no
test -n "$HDF5_LIBS" && hdf5_have=yes
nc_c_have=no
test -n "$NETCDF_C_LIBS" && nc_c_have=yes
nc_f_have=no
test -n "$NETCDFF_LIBS" && nc_f_have=yes
AC_MSG_NOTICE([HDF5: $hdf5_have])
AC_MSG_NOTICE([NetCDF C: $nc_c_have])
AC_MSG_NOTICE([NetCDF Fortran: $nc_f_have])
#AC_MSG_NOTICE([HDF5: $(test -n "$HDF5_LIBS" && echo yes || echo no)])
#AC_MSG_NOTICE([NetCDF C: $(test -n "$NETCDF_C_LIBS" && echo yes || echo no)])
#AC_MSG_NOTICE([NetCDF Fortran: $(test -n "$NETCDFF_LIBS" && echo yes || echo no)])
AC_MSG_NOTICE([Zoltan: $HAVE_ZOLTAN])
AC_MSG_NOTICE([VTK: $VTK at $with_vtk])
AC_MSG_NOTICE([Python embed: $HAVE_PYTHON])

dnl -------- Files to generate --------
AC_CONFIG_FILES([Makefile
          debug/Makefile
          bathymetry/Makefile
          ocean_forcing/Makefile ocean_forcing/tests/Makefile
          sediments/Makefile
          population_balance/Makefile
          hyperlight/Makefile
          femtools/Makefile femtools/tests/Makefile
          forward_interfaces/Makefile
          horizontal_adaptivity/Makefile horizontal_adaptivity/tests/Makefile
          preprocessor/Makefile
          error_measures/Makefile error_measures/tests/Makefile
          parameterisation/Makefile parameterisation/tests/Makefile
          fldecomp/Makefile
          assemble/Makefile assemble/tests/Makefile
          diagnostics/Makefile
          main/Makefile
          tools/Makefile
          tools/version-info
          python/setup.py
          climatology/Makefile
          libmba2d/Makefile
          libmba3d/Makefile
          libjudy/Makefile
          libjudy/src/Makefile
          libjudy/src/JudyCommon/Makefile
          libjudy/src/Judy1/Makefile
          libjudy/src/JudyL/Makefile
          libjudy/src/JudySL/Makefile
          libjudy/src/JudyHS/Makefile
          libwm/Makefile
          libvtkfortran/Makefile
          tests/tools.mk
          ICFERST/Makefile
          ICFERST/src/Makefile
          ICFERST/unittest/Makefile])

AC_OUTPUT
