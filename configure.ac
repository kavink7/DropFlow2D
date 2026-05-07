# ─────────────────────────────────────────────────────────────────────────────
# configure.ac  —  Autoconf input for dropflow2D
#
# Run via:  ./configure.sh   (bootstraps + configures in one step)
# Or manually:
#   autoreconf -fi
#   ./configure [options]
#
# Useful configure options
#   --prefix=<dir>      Install prefix  (default /usr/local)
#   --enable-debug      Build with -g -O0 instead of -O2
# ─────────────────────────────────────────────────────────────────────────────

# ── Package identity ─────────────────────────────────────────────────────────
AC_INIT([dropflow2D], [0.1.0], [], [dropflow2D])

# ── Automake options ─────────────────────────────────────────────────────────
# 'foreign'  : don't require NEWS/README/AUTHORS/ChangeLog
# 'subdir-objects': put .o files next to their .cpp (good hygiene)
AM_INIT_AUTOMAKE([foreign subdir-objects -Wall])

# ── Require a C++ compiler ───────────────────────────────────────────────────
AC_PROG_CXX
AC_LANG([C++])

# ── Enforce C++17 ────────────────────────────────────────────────────────────
# Try -std=c++17; fall back to -std=c++14 and error out if neither works.
# GCC >= 7 and Clang >= 5 support C++17 fully.
AX_CXX_STD=""
save_CXXFLAGS="$CXXFLAGS"
AC_MSG_CHECKING([whether $CXX accepts -std=c++17])
CXXFLAGS="$save_CXXFLAGS -std=c++17"
AC_COMPILE_IFELSE(
    [AC_LANG_PROGRAM([], [[]])],
    [AC_MSG_RESULT([yes]); AX_CXX_STD="-std=c++17"],
    [
        AC_MSG_RESULT([no])
        AC_MSG_CHECKING([whether $CXX accepts -std=c++14])
        CXXFLAGS="$save_CXXFLAGS -std=c++14"
        AC_COMPILE_IFELSE(
            [AC_LANG_PROGRAM([], [[]])],
            [AC_MSG_RESULT([yes]); AX_CXX_STD="-std=c++14"
             AC_MSG_WARN([C++17 not available — falling back to C++14.
             Some features may be unavailable. Upgrade to GCC >= 7 or Clang >= 5.])],
            [AC_MSG_RESULT([no])
             AC_MSG_ERROR([C++17 (or C++14) support is required but your compiler does not support it.
             Please upgrade your compiler.])]
        )
    ]
)
CXXFLAGS="$save_CXXFLAGS"   # restored; composed cleanly below

# ── Optional debug build ─────────────────────────────────────────────────────
AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug],
        [Build with -g -O0 instead of -O2 (default: no)])],
    [enable_debug=$enableval],
    [enable_debug=no])

AS_IF([test "x$enable_debug" = "xyes"],
    [OPT_FLAGS="-g -O0"; AC_MSG_NOTICE([Debug build enabled])],
    [OPT_FLAGS="-O2";    AC_MSG_NOTICE([Optimised build (-O2)  (use --enable-debug for -O0 -g)])])

# ── Compose final CXXFLAGS ───────────────────────────────────────────────────
# Preserve any user-supplied flags from the environment.
CXXFLAGS="$AX_CXX_STD $OPT_FLAGS -Wall -Wextra ${CXXFLAGS}"
AC_SUBST([CXXFLAGS])

# ── Check for required standard headers ──────────────────────────────────────
AC_CHECK_HEADERS([cmath vector string fstream sstream iostream cstdlib],
    [],
    [AC_MSG_WARN([A standard header was not found — build may fail])])

# ── sys/stat.h (POSIX, used for mkdir in output.cpp) ─────────────────────────
AC_CHECK_HEADER([sys/stat.h],
    [],
    [AC_MSG_ERROR([sys/stat.h is required but was not found])])

# ── Generate output files ─────────────────────────────────────────────────────
AC_CONFIG_FILES([Makefile])
AC_OUTPUT

# ── Summary ──────────────────────────────────────────────────────────────────
echo ""
echo "─────────────────────────────────────────────"
echo "  dropflow2D configuration summary"
echo ""
echo "  CXX      : $CXX"
echo "  C++ std  : $AX_CXX_STD"
echo "  CXXFLAGS : $CXXFLAGS"
echo "  prefix   : $prefix"
echo "  debug    : $enable_debug"
echo ""
echo "  Now run:  make"
echo "            make install   (optional)"
echo "─────────────────────────────────────────────"
