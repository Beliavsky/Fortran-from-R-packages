#!/usr/bin/env sh
set -eu
rm -rf build-gfortran
mkdir -p build-gfortran/mod build-gfortran/obj build-gfortran/bin
cd build-gfortran/obj
gfortran -std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface \
  -Werror=implicit-interface -fcheck=all -J../mod -I../mod -c \
  ../../src/quadprog_kinds.f90 ../../src/quadprog_core.f90 \
  ../../src/quadprog.f90 ../../src/bb_kinds.f90 ../../src/bb_interfaces.f90 \
  ../../src/bb_types.f90 ../../src/bb_projection.f90 ../../src/bb_spg.f90 \
  ../../src/bb_aux_optim.f90 ../../src/bb_nonlinear.f90 \
  ../../src/bb_drivers.f90 ../../src/bb.f90
ar rcs ../libbb.a ./*.o
cd ../..
for src in test/*.f90; do
  exe="build-gfortran/bin/$(basename "$src" .f90)"
  gfortran -std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface \
    -Werror=implicit-interface -fcheck=all -Ibuild-gfortran/mod \
    "$src" build-gfortran/libbb.a -o "$exe"
  "$exe"
done
