#!/usr/bin/env sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
BUILD="$ROOT/build-gfortran"
rm -rf "$BUILD"
mkdir -p "$BUILD"
cd "$BUILD"
FLAGS="-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0 -J."

gfortran $FLAGS -c "$ROOT/vendor/quadprog-fortran/src/quadprog_kinds.f90"
gfortran $FLAGS -c "$ROOT/vendor/quadprog-fortran/src/quadprog_core.f90"
gfortran $FLAGS -c "$ROOT/vendor/quadprog-fortran/src/quadprog.f90"
gfortran $FLAGS -c "$ROOT/src/quadprogxt.f90"
gfortran $FLAGS -c "$ROOT/test/qpxt_test_support.f90"

OBJS="quadprog_kinds.o quadprog_core.o quadprog.o quadprogxt.o qpxt_test_support.o"
for SRC in "$ROOT"/test/test_*.f90; do
  NAME=$(basename "$SRC" .f90)
  gfortran $FLAGS "$SRC" $OBJS -o "$NAME"
  "./$NAME"
done

for SRC in "$ROOT"/example/*.f90; do
  NAME=$(basename "$SRC" .f90)
  gfortran $FLAGS "$SRC" quadprog_kinds.o quadprog_core.o quadprog.o \
    quadprogxt.o -o "$NAME"
  "./$NAME"
done
