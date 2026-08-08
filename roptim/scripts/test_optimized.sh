#!/usr/bin/env sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
BUILD="$ROOT/build-gfortran-release"
rm -rf "$BUILD"
mkdir -p "$BUILD"
cd "$BUILD"
FLAGS="-std=f2018 -O3"
gfortran $FLAGS -c "$ROOT/src/roptim_lbfgsb_core.f90"
gfortran $FLAGS -c "$ROOT/src/roptim.f90"
gfortran $FLAGS -c "$ROOT/test/test_roptim.f90"
gfortran $FLAGS -o test_roptim roptim_lbfgsb_core.o roptim.o test_roptim.o
./test_roptim
gfortran $FLAGS -o rosenbrock_methods roptim_lbfgsb_core.o roptim.o \
  "$ROOT/example/rosenbrock_methods.f90"
./rosenbrock_methods
gfortran $FLAGS -o wild_sann roptim_lbfgsb_core.o roptim.o \
  "$ROOT/example/wild_sann.f90"
./wild_sann
