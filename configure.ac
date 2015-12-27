dnl $Id$
dnl
dnl
dnl Exported and configured variables:
dnl CFLAGS
dnl LDFLAGS
dnl INCLUDES


dnl ------------------------------------------------------------------
dnl Initialization for autoconf
dnl ------------------------------------------------------------------
AC_PREREQ(2.69)
AC_INIT([crfsuite], [0.13], [https://github.com/WladimirSidorenko])
AC_CONFIG_SRCDIR([frontend/main.c])
AC_CONFIG_MACRO_DIR([m4])

dnl ------------------------------------------------------------------
dnl Checks for system
dnl ------------------------------------------------------------------
AC_CANONICAL_HOST
AC_AIX
AC_MINIX
AC_ISC_POSIX


dnl ------------------------------------------------------------------
dnl Initialization for automake
dnl ------------------------------------------------------------------
AM_INIT_AUTOMAKE([subdir-objects parallel-tests color-tests])
AC_CONFIG_HEADERS(config.h)
AM_MAINTAINER_MODE
dnl AM_C_PROTOTYPES


dnl ------------------------------------------------------------------
dnl Checks for program
dnl ------------------------------------------------------------------
AM_PROG_CC_C_O
AC_PROG_CXX
AC_PROG_LIBTOOL
AC_PROG_INSTALL
AC_PROG_LN_S
dnl AX_PKG_SWIG
dnl AX_SWIG_PYTHON

dnl ------------------------------------------------------------------
dnl Initialization for variables
dnl ------------------------------------------------------------------
CFLAGS="-Wall -std=c99 ${ac_save_CFLAGS}"
LDFLAGS="${ac_save_LDFLAGS}"
dnl update `swig/python/setup.py.in` after changing this variable
CPPFLAGS="-I\${abs_top_srcdir} -I\${abs_top_srcdir}/include dnl
-I\${abs_top_srcdir}/lib/cqdb/include -I\${abs_srcdir}"
CXXFLAGS="-std=c++11"

dnl ------------------------------------------------------------------
dnl Checks for header files.
dnl ------------------------------------------------------------------
AC_HEADER_STDC
AC_CHECK_HEADERS(fcntl.h limits.h malloc.h strings.h unistd.h stdint.h)


dnl ------------------------------------------------------------------
dnl Checks for typedefs, structures, and compiler characteristics.
dnl ------------------------------------------------------------------
AC_C_CONST
AC_CHECK_SIZEOF
AC_TYPE_SIZE_T
AC_STRUCT_TM
AC_CHECK_SIZEOF(short)
AC_CHECK_SIZEOF(unsigned short)
AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(unsigned int)
AC_CHECK_SIZEOF(long)
AC_CHECK_SIZEOF(unsigned long)

AC_CHECK_TYPES([uint8_t, uint16_t, uint32_t])

dnl ------------------------------------------------------------------
dnl Checks for debugging mode
dnl ------------------------------------------------------------------
AC_ARG_ENABLE(
  debug,
  [AS_HELP_STRING([--enable-debug],[turn on debugging])]
)

if test "x$enable_debug" = "xyes"; then
   CFLAGS="-DDEBUG -O -g ${CFLAGS}"
else
   CFLAGS="-O3 -fomit-frame-pointer -ffast-math -Winline ${CFLAGS}"
fi

dnl ------------------------------------------------------------------
dnl Checks for profiling mode
dnl ------------------------------------------------------------------
AC_ARG_ENABLE(
  profile,
  [AS_HELP_STRING([--enable-profile],[turn on profiling])]
)

if test "x$enable_profile" = "xyes"; then
   CFLAGS="-DPROFILE -pg ${CFLAGS}"
fi


dnl ------------------------------------------------------------------
dnl Checks for SSE2 build
dnl ------------------------------------------------------------------
AC_ARG_ENABLE([sse2],
    AS_HELP_STRING(
        [--disable-sse2],
        [disable SSE2 optimization routines]
        )
    )

AS_IF([test "x$enable_sse2" != "xno"], [
    CFLAGS="-mfpmath=sse -msse2 -DUSE_SSE ${CFLAGS}"
])


dnl ------------------------------------------------------------------
dnl Checks for library functions.
dnl ------------------------------------------------------------------
AC_FUNC_ALLOCA
AC_FUNC_MEMCMP
AC_FUNC_VPRINTF
AC_CHECK_FUNCS(strdup strerror strtol strtoul)

dnl Check for math library
AC_CHECK_LIB(m, rand)

AC_ARG_WITH(
	liblbfgs,
	[AS_HELP_STRING([--with-liblbfgs=DIR],[liblbfgs directory])],
	[CPPFLAGS="${CPPFLAGS} -I${withval}/include"; LDFLAGS="${LDFLAGS} -L${withval}/lib"]
)
AC_CHECK_LIB(lbfgs, lbfgs)


dnl ------------------------------------------------------------------
dnl Set up test harness
dnl ------------------------------------------------------------------
AC_REQUIRE_AUX_FILE([tap-driver.sh])

dnl ------------------------------------------------------------------
dnl Export variables
dnl ------------------------------------------------------------------
AC_SUBST(CFLAGS)
AC_SUBST(LDFLAGS)
AC_SUBST(CPPFLAGS)
AC_SUBST(CXXFLAGS)
AC_SUBST(abs_top_srcdir)
AC_SUBST(includedir)
AC_SUBST(libdir)

dnl ------------------------------------------------------------------
dnl Output the configure results.
dnl ------------------------------------------------------------------
AC_CONFIG_FILES(Makefile genbinary.sh include/Makefile lib/cqdb/Makefile dnl
			 lib/crf/Makefile frontend/Makefile swig/Makefile dnl
			 tests/Makefile)
AC_CONFIG_FILES([swig/python/setup.py], [chmod +x swig/python/setup.py])

AC_OUTPUT
