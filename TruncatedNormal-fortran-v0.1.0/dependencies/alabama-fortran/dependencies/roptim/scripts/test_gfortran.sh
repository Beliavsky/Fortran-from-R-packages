#!/usr/bin/env sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
BUILD="$ROOT/build-gfortran-check"
rm -rf "$BUILD"
mkdir -p "$BUILD"
cd "$BUILD"
FLAGS="-std=f2018 -O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -Wimplicit-interface -Werror=implicit-interface"
gfortran $FLAGS -c "$ROOT/src/roptim_lbfgsb_core.f90"
gfortran $FLAGS -c "$ROOT/src/roptim.f90"
gfortran $FLAGS -c "$ROOT/test/test_roptim.f90"
gfortran $FLAGS -o test_roptim roptim_lbfgsb_core.o roptim.o test_roptim.o
./test_roptim
