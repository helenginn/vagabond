#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(vagabond, 0.0, vagabond@hginn.co.uk)
AM_INIT_AUTOMAKE()
LT_INIT()

#AC_CONFIG_SRCDIR([])
#AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
#AC_PROG_CC

# Checks for libraries.
AC_CHECK_LIB([png],[png_read_info],[],[AC_MSG_ERROR([Cannot find libpng])])

# Checks for header files.
AC_CHECK_HEADERS([fcntl.h float.h limits.h malloc.h stddef.h stdlib.h string.h sys/file.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_MODE_T
AC_TYPE_OFF_T
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ALLOCA
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor getcwd memmove memset mkdir pow putenv rint setenv sqrt strchr strdup strerror strrchr strstr])

AC_CONFIG_SUBDIRS([libfftw])

AC_OUTPUT(Makefile libinfo/Makefile libccp4/Makefile libica/Makefile libsrc/Makefile) #libfftw/Makefile)

