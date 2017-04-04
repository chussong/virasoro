AC_INIT(virasoro, 1.0.2)
AM_INIT_AUTOMAKE([foreign])

AC_ARG_ENABLE([dev],
			  AS_HELP_STRING([--enable-dev], [Enable additional warnings and errors.]))

AC_ARG_ENABLE([wstp],
			  AS_HELP_STRING([--disable-wstp], [Do not build support for calling from Mathematica.]))

AC_ARG_ENABLE([verbose],
			  AS_HELP_STRING([--enable-verbose], [Enable display of intermediate results.]))

dnl Switch to a C++ compiler, and check if it works.
AC_LANG(C++)
dnl ${CFLAGS="-O3"}
AC_PROG_CC
dnl ${CXXFLAGS="-O3"}
AC_PROG_CXX
AC_CANONICAL_HOST

dnl Check for GMP and derivatives.
AC_CHECK_HEADERS(gmp.h, [], [AC_ERROR(["GMP not detected. Either you have not installed it or the directory where it's installed is not set up to be searched."])])
dnl gmplib=/usr/local/lib/libgmp.a
AC_CHECK_HEADERS(mpfr.h, [], [AC_ERROR(["MPFR not detected. Either you have not installed it or the directory where it's installed is not set up to be searched."])])
dnl mpfrlib=/usr/local/lib/libmpfr.a
AC_CHECK_HEADERS(mpc.h, [], [AC_ERROR(["MPC not detected. Either you have not installed it or the directory where it's installed is not set up to be searched."])])
dnl mpclib=/usr/local/lib/libmpc.a

dnl headers should probably be split for the different targets
headers="mpreal.h mpcomplex.h cpqmn.h hmn.h runfile.h virasoro.h vwstp.h access.h"
object_files="virasoro.o runfile.o"
cflags="-O3"
cxxflags="-O3 -std=c++14"
devflags="-Wall -Wextra -pedantic -Werror -g"
ldadd="-lmpc -lmpfr -lgmp -lpthread"
dnl staticlibs="$mpclib $mpfrlib $gmplib"

dnl Check for WSTP
build_linux=no
build_mac=no

case "${host_os}" in
	linux*)
		AC_MSG_NOTICE([Building for Linux.])
		build_linux=yes
		if test "${MATHEMATICA_HOME+set}" = set; then
			AC_MSG_NOTICE([Mathematica found via environment variable.])
			wstpdir=$MATHEMATICA_HOME/SystemFiles/Links/WSTP/DeveloperKit/Linux-x86-64/CompilerAdditions
			AC_CHECK_FILE($wstpdir/libWSTP64i4.a,
						  [wstp=1],
						  [wstpdir=$MATHEMATICA_HOME/SystemFiles/Links/WSTP/DeveloperKit/Linux/CompilerAdditions
						   AC_CHECK_FILE($wstpdir/libWSTP64i4.a,
										 [wstp=1],
										 [wstp=0])])
			if test $wstp == 1
			then
				wstplib=$wstpdir/libWSTP64i4.a
				wstpextralibs="$wstplib -lm -lrt -lstdc++ -ldl -luuid"
				AC_SUBST([pkgdatadir], [$HOME/.Mathematica/Applications])
			fi
		else
			wstp=0
		fi
		;;
	darwin*)
		AC_MSG_NOTICE([Building for Mac OS.])
		build_mac=yes
		AC_CHECK_FILE(/Applications/Mathematica.app/Contents/SystemFiles/Links/WSTP/DeveloperKit/MacOSX-x86-64/CompilerAdditions/wsprep, [wstp=1], [wstp=0])
		if test $wstp == 1
		then
			AC_MSG_NOTICE([Mathematica found in the expected location.])
			wstpdir=/Applications/Mathematica.app/Contents/SystemFiles/Links/WSTP/DeveloperKit/MacOSX-x86-64/CompilerAdditions
			wstplib=$wstpdir/libWSTPi4.a
			wstpextralibs="$wstplib -lc++ -framework Foundation"
			AC_SUBST([pkgdatadir], [$HOME/Library/Mathematica/Applications])
		fi
		;;
	*)
		AC_MSG_NOTICE([Usable build OS not detected.])
		;;
esac

AS_IF([test "x$enable_wstp" = "xno"], [
	AC_MSG_NOTICE([WSTP has been manually disabled.])
	wstp=0
])

if test $wstp == 1; then
	AC_MSG_NOTICE([Telling make to build VWSTP.])
	VWSTP="vwstp"
	AC_DEFINE([HAVE_WSTP], [1], "Enable WSTP functions inside shared sources.")
	AC_SUBST([VWSTP], ['vwstp${EXEEXT}'])
	wstpflags="-Wno-unused-parameter -I$wstpdir"
	object_files="$object_files vwstp.o"
	wsprep=$wstpdir/wsprep
	AC_SUBST([USERNAME], [$USER])
else
	VWSTP=""
	object_files="$object_files main.o"
fi
AM_CONDITIONAL([MAKE_VWSTP], [test $wstp = 1])

AS_IF([test "x$enable_dev" = "xyes"], [
	cflags="$cflags $devflags"
	cxxflags="$cxxflags $devflags"
])

AS_IF([test "x$enable_verbose" = "xyes"], [
	cxxflags="$cxxflags -DVERBOSE_DEBUG"
])

dnl MAKE SURE CFLAGS, LDFLAGS, STATICLIBS GET SENT TO THE MAKEFILE AND ARE USED IN CORRECT PROGS
AC_SUBST([AM_CFLAGS], [$cflags])
AC_SUBST([AM_CXXFLAGS], [$cxxflags])
AC_SUBST([AM_LDFLAGS], [$ldflags])
AC_SUBST([AM_LDADD], [$ldadd])
AC_SUBST([WSPREP],[$wsprep])
AC_SUBST([VWSTP])
AC_SUBST([wstpflags])
AC_SUBST([wstpextralibs])
dnl AC_SUBST([staticlibs])
AC_SUBST([object_files])
AC_CONFIG_HEADERS([include/config.h])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
