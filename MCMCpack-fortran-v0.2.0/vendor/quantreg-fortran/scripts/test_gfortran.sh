#!/usr/bin/env sh
set -eu
FC=${FC:-gfortran}
FLAGS='-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0'
rm -rf build-strict
mkdir build-strict
$FC $FLAGS -J build-strict -I build-strict -c \
  src/quantreg_kinds.f90 src/quantreg_types.f90 src/quantreg_linalg.f90 \
  src/quantreg_dense.f90 src/quantreg_select.f90 src/quantreg_local.f90 \
  src/quantreg_nonlinear.f90 src/quantreg_bootstrap.f90 \
  src/quantreg_utils.f90 src/quantreg.f90
OBJS='quantreg_kinds.o quantreg_types.o quantreg_linalg.o quantreg_dense.o quantreg_select.o quantreg_local.o quantreg_nonlinear.o quantreg_bootstrap.o quantreg_utils.o quantreg.o'
for t in test/*.f90; do
  e=build-strict/$(basename "$t" .f90)
  $FC $FLAGS -I build-strict "$t" $OBJS -o "$e"
  "$e"
done
for t in example/*.f90; do
  e=build-strict/$(basename "$t" .f90)
  $FC $FLAGS -I build-strict "$t" $OBJS -o "$e"
  "$e"
done
