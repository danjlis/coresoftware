AC_INIT(pgcal,[1.00])
AC_CONFIG_SRCDIR([configure.ac])

AM_INIT_AUTOMAKE

dnl order is important, if AC_PROG_LIBTOOL is called before AC_PROG_CXX(CC g++)
dnl it will choose the g++ compiler over CC on Solaris
AC_PROG_CXX(CC g++ clang++)
LT_INIT([disable-static])

dnl   no point in suppressing warnings people should
dnl   at least see them, so here we go for g++: -Wall
CXXFLAGS="$CXXFLAGS -Wall -Wextra -Wshadow -Werror"

CINTDEFS=" -noIncludePaths  -inlineInputHeader "
AC_SUBST(CINTDEFS)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
