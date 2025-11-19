AC_DEFUN([CHECK_ZOLTAN], [
    AS_IF([test -n "$with_zoltan"], [
      CPPFLAGS="$CPPFLAGS -I$with_zoltan/include"
      FCFLAGS="$FCFLAGS -I$with_zoltan/include"
      LDFLAGS="$LDFLAGS -L$with_zoltan/lib"
    ])
    AC_CHECK_HEADERS([zoltan.h],[],[
    have_zoltan_hdr=no
    # Try with Trilinos include path if available
    unset ac_cv_header_zoltan_h
    save_CPPFLAGS="$CPPFLAGS"
    CPPFLAGS="$CPPFLAGS -I/usr/include/trilinos"
    AC_CHECK_HEADERS([zoltan.h], [have_zoltan_hdr=yes], [
      CPPFLAGS="$save_CPPFLAGS"
    ])
  ])
  AC_MSG_CHECKING([whether to link against zoltan or trilinos_zoltan])
  have_zoltan_lib=no
  AC_CHECK_LIB([trilinos_zoltan],[Zoltan_Create],[have_zoltan_lib=yes; ZOLTAN_LIBS="-ltrilinos_zoltan"],[
    AC_CHECK_LIB([zoltan],[Zoltan_Create],[have_zoltan_lib=yes; ZOLTAN_LIBS="-lzoltan"],[
      AC_MSG_WARN([Neither zoltan nor trilinos_zoltan library found])
    ])
  ])
  AC_SUBST([ZOLTAN_LIBS])
  LIBS="$LIBS $ZOLTAN_LIBS"

  AS_IF([test "x$have_zoltan_hdr" != "xno" -a "x$have_zoltan_lib" != "xno"],[
    AC_LANG_PUSH([Fortran])
    AC_MSG_CHECKING([for usable zoltan Fortran module])
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM([], [
      use zoltan
      integer :: i
      i = 0
  ])],
  [AC_DEFINE([HAVE_ZOLTAN],[1],[Have Zoltan])
   AC_DEFINE([HAVE_ZOLTAN_MOD],[1],[Have Zoltan Fortran module])
   HAVE_ZOLTAN=yes
   AC_MSG_RESULT([yes])
],
  [
    save_FCFLAGS="$FCFLAGS"
    FCFLAGS="$FCFLAGS -I/usr/include"
    AC_COMPILE_IFELSE(
      [AC_LANG_PROGRAM([], [
          use zoltan
          integer :: i
          i = 0
      ])],
      [AC_DEFINE([HAVE_ZOLTAN],[1],[Have Zoltan])
       AC_DEFINE([HAVE_ZOLTAN_MOD],[1],[Have Zoltan Fortran module])
       HAVE_ZOLTAN=yes
       AC_MSG_RESULT([yes])
      ],
      [
        AC_MSG_RESULT([no])
        FCFLAGS="$save_FCFLAGS"
      ]
    )
]
)
    AC_LANG_POP([Fortran])
    ],[:])
])
