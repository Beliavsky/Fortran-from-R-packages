#!/bin/bash
#
# build_scip.sh — Build SoPlex + SCIP static libraries from source
#
# Called by configure during R CMD INSTALL.
# Produces: src/sciplib/ with include/ and lib/ subdirectories.
#
# Modeled on inst/build_highs.sh from the highs R package.

#
# Detect tools
#
if test -z "${MAKE}"; then MAKE=`which make 2>/dev/null`; fi
if test -z "${MAKE}"; then MAKE=`which /Applications/Xcode.app/Contents/Developer/usr/bin/make 2>/dev/null`; fi

if test -z "${CMAKE_EXE}"; then CMAKE_EXE=`which cmake4 2>/dev/null`; fi
if test -z "${CMAKE_EXE}"; then CMAKE_EXE=`which cmake3 2>/dev/null`; fi
if test -z "${CMAKE_EXE}"; then CMAKE_EXE=`which cmake 2>/dev/null`; fi
if test -z "${CMAKE_EXE}"; then CMAKE_EXE=`which /Applications/CMake.app/Contents/bin/cmake 2>/dev/null`; fi

if test -z "${CMAKE_EXE}"; then
    echo "Could not find 'cmake'!"
    exit 1
fi

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
    echo "'R_HOME' could not be found!"
    exit 1
fi

#
# Get compiler settings from R
#
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`

export CC=`"${R_HOME}/bin/R" CMD config CC`
export CXX=`"${R_HOME}/bin/R" CMD config CXX17`

R_SCIP_PKG_HOME=`pwd`
CONFIG_DIR=${R_SCIP_PKG_HOME}/inst/config

# r_streams.h lives in CONFIG_DIR — SCIP/SoPlex r_pkg branches include it
# SCIP_LONGINT_FORMAT: override SCIP's default "I64d" on Windows (MSVC format)
# to "lld" which works with MinGW/gcc in Rtools
export CFLAGS="${CFLAGS} -I${CONFIG_DIR} -DSCIP_LONGINT_FORMAT=\\\"lld\\\""
export CXXFLAGS="${CXXFLAGS} -I${CONFIG_DIR} -DSCIP_LONGINT_FORMAT=\\\"lld\\\""
export LDFLAGS

echo ""
echo "CMAKE VERSION: '`${CMAKE_EXE} --version | head -n 1`'"
echo "CC: '${CC}'"
echo "CXX: '${CXX}'"
echo "CFLAGS: '${CFLAGS}'"
echo "CXXFLAGS: '${CXXFLAGS}'"
echo ""

SOPLEX_SRC_DIR=${R_SCIP_PKG_HOME}/inst/soplex
SCIP_SRC_DIR=${R_SCIP_PKG_HOME}/inst/scip
SOPLEX_INSTALL_DIR=${R_SCIP_PKG_HOME}/src/soplexlib
SCIP_INSTALL_DIR=${R_SCIP_PKG_HOME}/src/sciplib

#
# Detect ccache for faster rebuilds
#
CCACHE_OPTS=""
CCACHE_EXE=`which ccache 2>/dev/null`
if test -n "${CCACHE_EXE}"; then
    CCACHE_OPTS="-DCMAKE_C_COMPILER_LAUNCHER=${CCACHE_EXE} -DCMAKE_CXX_COMPILER_LAUNCHER=${CCACHE_EXE}"
    echo "Found ccache: ${CCACHE_EXE}"
fi

#
# Common CMake options
#
COMMON_CMAKE_OPTS="
    -DCMAKE_POSITION_INDEPENDENT_CODE:bool=ON
    -DBUILD_SHARED_LIBS:bool=OFF
    -DBUILD_TESTING:bool=OFF
    -DCMAKE_VERBOSE_MAKEFILE:bool=ON
    ${CCACHE_OPTS}
"

# Platform-specific CMake flag
if test "$(uname -s)" = "Darwin"; then
    CMAKE_PLATFORM_OPTS="-DCMAKE_HOST_APPLE:bool=ON"
else
    CMAKE_PLATFORM_OPTS="-G \"Unix Makefiles\""
fi

# ========================================================
# Step 1: Build SoPlex (static library only)
# ========================================================
echo ">>> Building SoPlex..."
SOPLEX_BUILD_DIR=${SOPLEX_SRC_DIR}/build
mkdir -p ${SOPLEX_BUILD_DIR}
mkdir -p ${SOPLEX_INSTALL_DIR}/lib
mkdir -p ${SOPLEX_INSTALL_DIR}/include/soplex
cd ${SOPLEX_BUILD_DIR}

# R CMD build strips directories named "check" — recreate stubs if missing
# (SoPlex CMakeLists.txt unconditionally does add_subdirectory(check))
if test ! -d ${SOPLEX_SRC_DIR}/check; then
    mkdir -p ${SOPLEX_SRC_DIR}/check
    echo "# stub" > ${SOPLEX_SRC_DIR}/check/CMakeLists.txt
fi

SOPLEX_CMAKE_OPTS="
    ${COMMON_CMAKE_OPTS}
    -DCMAKE_INSTALL_PREFIX=${SOPLEX_INSTALL_DIR}
    -DZLIB:bool=OFF
    -DGMP:bool=OFF
    -DMPFR:bool=OFF
    -DBOOST:bool=OFF
    -DPAPILO:bool=OFF
    -DQUADMATH:bool=OFF
"

eval ${CMAKE_EXE} .. ${SOPLEX_CMAKE_OPTS} ${CMAKE_PLATFORM_OPTS} || exit 1

# Build only the static library target (skip shared lib which fails
# due to unresolved R symbols like Rprintf, r_cout, r_cerr)
${MAKE} libsoplex || exit 1

echo ">>> SoPlex built in ${SOPLEX_BUILD_DIR}"

# Manual install: copy static lib + headers
cp lib/libsoplex.a ${SOPLEX_INSTALL_DIR}/lib/
cp -r ${SOPLEX_SRC_DIR}/src/soplex/*.h ${SOPLEX_INSTALL_DIR}/include/soplex/
cp soplex/config.h ${SOPLEX_INSTALL_DIR}/include/soplex/

echo ">>> SoPlex installed to ${SOPLEX_INSTALL_DIR}"

# ========================================================
# Step 2: Build SCIP
# ========================================================
echo ">>> Building SCIP..."
SCIP_BUILD_DIR=${SCIP_SRC_DIR}/build
mkdir -p ${SCIP_BUILD_DIR}
mkdir -p ${SCIP_INSTALL_DIR}/lib
mkdir -p ${SCIP_INSTALL_DIR}/include
cd ${SCIP_BUILD_DIR}

# SOPLEX_DIR points at the SoPlex *build* dir where soplex-config.cmake lives
SCIP_CMAKE_OPTS="
    ${COMMON_CMAKE_OPTS}
    -DCMAKE_INSTALL_PREFIX=${SCIP_INSTALL_DIR}
    -DSOPLEX_DIR=${SOPLEX_BUILD_DIR}
    -DZLIB:bool=OFF
    -DGMP:bool=OFF
    -DREADLINE:bool=OFF
    -DPAPILO:bool=OFF
    -DZIMPL:bool=OFF
    -DAMPL:bool=OFF
    -DIPOPT:bool=OFF
    -DWORHP:bool=OFF
    -DCONOPT:bool=OFF
    -DLAPACK:bool=OFF
    -DAUTOBUILD:bool=OFF
    -DSHARED:bool=OFF
    -DEXPRINT=none
    -DLPS=spx
    -DTPI=${TPI_MODE:-none}
"

eval ${CMAKE_EXE} .. ${SCIP_CMAKE_OPTS} ${CMAKE_PLATFORM_OPTS} || exit 1

# Build only the static library target
${MAKE} libscip || exit 1

# Manual install: copy static lib + all headers preserving directory structure
# SCIP headers use paths like <scip/scip.h>, <blockmemshell/memory.h>, <lpi/type_lpi.h>
cp lib/libscip.a ${SCIP_INSTALL_DIR}/lib/
# Copy the entire src/ tree (headers only) preserving subdirectory structure
cd ${SCIP_SRC_DIR}/src
find . -name '*.h' | while read f; do
    dir="${SCIP_INSTALL_DIR}/include/$(dirname "$f")"
    mkdir -p "$dir"
    cp "$f" "$dir/"
done
cd ${SCIP_BUILD_DIR}
# Copy CMake-generated headers (config.h, scip_export.h)
mkdir -p ${SCIP_INSTALL_DIR}/include/scip
cp scip/config.h ${SCIP_INSTALL_DIR}/include/scip/
cp scip/scip_export.h ${SCIP_INSTALL_DIR}/include/scip/

echo ">>> SCIP installed to ${SCIP_INSTALL_DIR}"

cd ${R_SCIP_PKG_HOME}
