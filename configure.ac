#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([bpopsim], [0.1], [jeffrey.e.barrick@gmail.com], [bpopsim])
AC_CONFIG_AUX_DIR(aux_build)
AC_CONFIG_MACRO_DIR([aux_build/m4])
AC_CONFIG_HEADERS([aux_build/config.h])
#AC_CONFIG_SRCDIR([])

LT_INIT
AM_INIT_AUTOMAKE

# Checks for programs.
AM_PROG_CC_C_O
AC_PROG_CXX
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_MAKE_SET

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([stdint.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([pow])

#checks for libraries                
        
## gsl
AX_PATH_GSL(1.14, [], [AC_MSG_ERROR(["libgsl is required."])] )
##AC_SEARCH_LIBS([gsl_rng_alloc], [gsl], [], [AC_MSG_ERROR(["libgsl >=1.14 is required."])])
##AC_CHECK_LIB([gsl], [gsl_rng_alloc])

## boost
#BOOST_REQUIRE([1.42])
#BOOST_PROGRAM_OPTIONS([s])

AC_SUBST([BPOPSIM], [src/c/bpopsim])

AC_CONFIG_FILES([
	src/c/bpopsim/Makefile
	Makefile
])
AC_OUTPUT


