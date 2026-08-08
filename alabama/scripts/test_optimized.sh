#!/usr/bin/env sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$root"
rm -rf build_optimized
mkdir -p build_optimized/mod
srcs="dependencies/numDeriv-fortran/src/numderiv_kinds.f90 dependencies/numDeriv-fortran/src/numderiv_types.f90 dependencies/numDeriv-fortran/src/numderiv_callbacks.f90 dependencies/numDeriv-fortran/src/numderiv_core.f90 dependencies/numDeriv-fortran/src/numderiv.f90 dependencies/roptim/src/roptim_lbfgsb_core.f90 dependencies/roptim/src/roptim.f90 src/alabama.f90"
gfortran -std=f2018 -O3 -DNDEBUG -Wimplicit-interface -Werror=implicit-interface -J build_optimized/mod -I build_optimized/mod $srcs test/test_alabama.f90 -o build_optimized/test_alabama
./build_optimized/test_alabama
