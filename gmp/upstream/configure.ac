dnl Process this file with autoconf to produce a configure script.
AC_INIT(src/biginteger.cc)

AC_ARG_WITH([gmp-include],
            AS_HELP_STRING([--with-gmp-include=INCLUDE_PATH],
                           [the location of GMP header files]),
            [gmp_include_path=$withval])

if test -n "$gmp_include_path" ; then
   GMP_CPPFLAGS="-I${gmp_include_path}"
elif test -n "${gmp_INCLUDE}" ; then
   GMP_CPPFLAGS="-I${gmp_INCLUDE}"
fi

dnl Ditto for a library path  environment variable  >> gmp_LDFLAGS <<
AC_ARG_WITH([gmp-lib],
            AS_HELP_STRING([--with-gmp-lib=LIB_PATH],
                           [the location of GMP libraries]),
            [gmp_lib_path=$withval])

if test -n "$gmp_lib_path" ; then
   GMP_LDFLAGS="-L$gmp_lib_path ${LDFLAGS}"
elif test -n "${gmp_LDFLAGS}" ; then
   GMP_LDFLAGS="-L${gmp_LDFLAGS} ${LDFLAGS}"
fi


CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
CXX=`"${R_HOME}/bin/R" CMD config CXX`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`
MAKE=`"${R_HOME}/bin/R" CMD config MAKE`
AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX


dnl use the include and lib
CPPFLAGS="${CPPFLAGS} ${MPFR_CPPFLAGS}"
LDFLAGS="${LDFLAGS} ${MPFR_LDFLAGS}"

PKG_CPPFLAGS=$CPPFLAGS
PKG_CFLAGS=$CFLAGS
PKG_LDFLAGS=$LDFLAGS



dnl Check for lib gmp

AC_CHECK_HEADER([gmp.h], ,
    [AC_MSG_ERROR([Header file gmp.h not found; maybe use --with-gmp-include=INCLUDE_PATH])])
AC_CHECK_LIB(gmp,__gmpz_ui_sub,, [AC_MSG_ERROR([GNU MP not found, or not 4.1.4 or up, see http://gmplib.org])])

AC_SUBST(PKG_CPPFLAGS)
AC_SUBST(PKG_CFLAGS)
AC_SUBST(PKG_LDFLAGS)

AC_CONFIG_FILES(src/Makevars)
AC_OUTPUT





