#!/usr/bin/env sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$root"
rm -rf build_gfortran
mkdir -p build_gfortran/mod
srcs="dependencies/numDeriv-fortran/src/numderiv_kinds.f90 dependencies/numDeriv-fortran/src/numderiv_types.f90 dependencies/numDeriv-fortran/src/numderiv_callbacks.f90 dependencies/numDeriv-fortran/src/numderiv_core.f90 dependencies/numDeriv-fortran/src/numderiv.f90 dependencies/roptim/src/roptim_lbfgsb_core.f90 dependencies/roptim/src/roptim.f90 src/alabama.f90"
gfortran -std=f2018 -O0 -g -Wimplicit-interface -Werror=implicit-interface -fcheck=all -ffpe-trap=invalid,zero,overflow -J build_gfortran/mod -I build_gfortran/mod $srcs test/test_alabama.f90 -o build_gfortran/test_alabama
./build_gfortran/test_alabama
gfortran -std=f2018 -O0 -g -Wimplicit-interface -Werror=implicit-interface -fcheck=all -ffpe-trap=invalid,zero,overflow -J build_gfortran/mod -I build_gfortran/mod $srcs example/constrained_example.f90 -o build_gfortran/constrained_example
./build_gfortran/constrained_example
