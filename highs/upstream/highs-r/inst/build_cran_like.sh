#!/bin/bash

# Script to reproduce CRAN-like build environment
# Based on the build log from CRAN's clang check

set -e

echo "========================================="
echo "CRAN-like Build Environment Setup"
echo "========================================="

# Set compilers to clang (adjust paths if needed)
export CC="clang -std=gnu23"
export CXX="clang++"
export CXX11="clang++"

# Set compiler flags matching CRAN's configuration
export CXXFLAGS="-O3 -Wall -pedantic -frtti -Wp,-D_FORTIFY_SOURCE=3 -Wno-missing-template-arg-list-after-template-kw"
export CFLAGS="-O3 -Wall -pedantic -Wp,-D_FORTIFY_SOURCE=3"

# Additional flags that may be needed
export CPPFLAGS=""
export LDFLAGS=""

echo "CC: ${CC}"
echo "CXX: ${CXX}"
echo "CXXFLAGS: ${CXXFLAGS}"
echo "CFLAGS: ${CFLAGS}"
echo "========================================="
echo ""

# Clean previous builds
echo "Cleaning previous builds..."
rm -f highs_*.tar.gz
rm -rf highs.Rcheck
rm -rf src/*.o src/*.so

# Build the package
echo "Building package..."
R CMD build .

# Optionally run check
if [ "$1" = "--check" ]; then
    echo ""
    echo "Running R CMD check --as-cran..."
    R CMD check --as-cran highs_*.tar.gz
fi

echo ""
echo "Build complete!"
echo "To install: R CMD INSTALL highs_*.tar.gz"
