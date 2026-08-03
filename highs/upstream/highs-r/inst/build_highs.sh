#!/bin/bash


#
# Set Variables
#
if test -z "${MAKE}"; then MAKE=`which make` 2> /dev/null; fi
if test -z "${MAKE}"; then MAKE=`which /Applications/Xcode.app/Contents/Developer/usr/bin/make` 2> /dev/null; fi


if test -z "${CMAKE_EXE}"; then CMAKE_EXE=`which cmake4` 2> /dev/null; fi
if test -z "${CMAKE_EXE}"; then CMAKE_EXE=`which cmake3` 2> /dev/null; fi
if test -z "${CMAKE_EXE}"; then CMAKE_EXE=`which cmake2` 2> /dev/null; fi
if test -z "${CMAKE_EXE}"; then CMAKE_EXE=`which cmake` 2> /dev/null; fi
if test -z "${CMAKE_EXE}"; then CMAKE_EXE=`which /Applications/CMake.app/Contents/bin/cmake` 2> /dev/null; fi

if test -z "${CMAKE_EXE}"; then
    echo "Could not find 'cmake'!"
    exit 1
fi


if test -z "${CCACHE}"; then CCACHE=`which ccache` 2> /dev/null; fi


: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
    echo "'R_HOME' could not be found!"
    exit 1
fi


# PKG_CPPFLAGS=`Rscript -e 'Rcpp:::CxxFlags()'`
# PKG_LIBS=`Rscript -e 'Rcpp:::LdFlags()'` 
# CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config --cppflags`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
R_INCLUDE_DIR=`"${R_HOME}/bin/Rscript" -e 'writeLines(R.home("include"))'`
R_INCLUDE_DIR="-I${R_INCLUDE_DIR}"
RCPP_INCLUDE_DIR=`"${R_HOME}/bin/Rscript" -e 'cat(system.file("include", package = "Rcpp"))'`
RCPP_INCLUDE_DIR="-I${RCPP_INCLUDE_DIR}"

export CC=`"${R_HOME}/bin/R" CMD config CC`
export CXX=`"${R_HOME}/bin/R" CMD config CXX17`
export CXX17=`"${R_HOME}/bin/R" CMD config CXX17`
export CPPFLAGS="${CPPFLAGS} ${RCPP_INCLUDE_DIR} ${R_INCLUDE_DIR}"
export CFLAGS="${CFLAGS} ${R_INCLUDE_DIR}"
export CXXFLAGS="${CXXFLAGS} ${RCPP_INCLUDE_DIR} ${R_INCLUDE_DIR}"
export LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`
export CXX_STD=`"${R_HOME}/bin/R" CMD config CXX_STD`

echo ""
echo "CMAKE VERSION: '`${CMAKE_EXE} --version | head -n 1`'"
echo "arch: '$(arch)'"
echo "CC: '${CC}'"
echo "CXX: '${CXX}'"
echo "CXX17: '${CXX}'"
echo "CXXFLAGS: '${CXXFLAGS}'"
echo "CFLAGS: '${CFLAGS}'"
echo "CPPFLAGS: '${CPPFLAGS}'"
echo "LDFLAGS: '${LDFLAGS}'"
echo "R_HIGHS_BUILD_DIR: '${R_HIGHS_BUILD_DIR}'"
echo "R_HIGHS_LIB_DIR: '${R_HIGHS_LIB_DIR}'"
echo "R_INCLUDE_DIR: '${R_INCLUDE_DIR}'"
echo ""


R_HIGHS_PKG_HOME=`pwd`
HIGHS_SRC_DIR=${R_HIGHS_PKG_HOME}/inst/HiGHS
R_HIGHS_BUILD_DIR=${HIGHS_SRC_DIR}/build
R_HIGHS_LIB_DIR=${R_HIGHS_PKG_HOME}/src/highslib


#
# Setup Build Directory
#
mkdir -p ${R_HIGHS_BUILD_DIR}
mkdir -p ${R_HIGHS_LIB_DIR}
cd ${R_HIGHS_BUILD_DIR}


#
# Derive Build Options
#

# The flag CMAKE_CXX_COMPILER_WORKS signals cmake that we know
# that the compiler works therefore cmake has not to test.
# Documentation on this flag is very hard to find.
# https://www.linuxfixes.com/2021/10/solved-cmake-c-compiler-is-not-able-to.html?m=0
# https://github.com/Kitware/CMake/blob/master/Modules/CMakeTestCXXCompiler.cmake

# In R-oldrelease HiGHS tries to build with 'NMake Makefiles' instead of 'Unix Makefiles'
# the additional flag '-G "Unix Makefiles"' forces the use of 'Unix Makefiles'.
# But than it fails since prior to 4.2.0 R-win tries to compile and test for 'i386' and 'x64'
# and fails for 'i386'.
# if test "${OS_TYPE}" = "unix"; then
DEFAULT_CMAKE_OPTS="
    -DCMAKE_INSTALL_PREFIX=${R_HIGHS_LIB_DIR} \
    -DCMAKE_POSITION_INDEPENDENT_CODE:bool=ON \
    -DBUILD_SHARED_LIBS:bool=OFF \
    -DBUILD_TESTING:bool=OFF \
    -DZLIB:bool=OFF \
    -DBUILD_EXAMPLES:bool=OFF \
    -DFAST_BUILD:bool=ON \
    -DCMAKE_VERBOSE_MAKEFILE:bool=ON \
    -DCUPDLP_GPU:bool=OFF \
    -DUSE_DOTNET_STD_21:bool=OFF \
"
COMPILER_LAUNCHER="-DCMAKE_C_COMPILER_LAUNCHER=${CCACHE} -DCMAKE_CXX_COMPILER_LAUNCHER=${CCACHE}"
if test -z "${CCACHE}"; then
    CMAKE_OPTS=${DEFAULT_CMAKE_OPTS}
else
    CMAKE_OPTS="${COMPILER_LAUNCHER} ${DEFAULT_CMAKE_OPTS}"
fi

echo ""
echo "Cmake Options"
echo ${CMAKE_OPTS}
echo ""


#
# Build and Install HiGHS
#
# Currently there would be no distinction necessary
if test "$(uname -s)" = "Darwin"; then
    ${CMAKE_EXE} .. ${CMAKE_OPTS} -DCMAKE_HOST_APPLE:bool=ON
else
    ${CMAKE_EXE} .. ${CMAKE_OPTS} -G "Unix Makefiles"
fi

${MAKE} install
