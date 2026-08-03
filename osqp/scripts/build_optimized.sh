#!/usr/bin/env sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
BUILD="$ROOT/build/optimized"
rm -rf "$BUILD"
mkdir -p "$BUILD"
cd "$BUILD"
CFLAGS='-std=c99 -Wall -Wextra -pedantic -O3'
FFLAGS='-std=f2018 -Wall -Wextra -Wimplicit-interface -O3'
gcc $CFLAGS -I"$ROOT/src" -c "$ROOT/src/osqp_dynamic_loader.c"
for src in osqp_kinds osqp_constants osqp_sparse osqp_c_bindings osqp_types osqp_solver_api osqp_model_api osqp; do
  gfortran $FFLAGS -J. -I. -c "$ROOT/src/$src.f90"
done
OBJS='osqp_dynamic_loader.o osqp_kinds.o osqp_constants.o osqp_sparse.o osqp_c_bindings.o osqp_types.o osqp_solver_api.o osqp_model_api.o osqp.o'
for src in "$ROOT"/test/*.f90; do
  exe=$(basename "$src" .f90)
  gfortran $FFLAGS -I. -o "$exe" $OBJS "$src" -ldl
  (cd "$ROOT" && "$BUILD/$exe")
done
for src in "$ROOT"/example/*.f90 "$ROOT"/app/*.f90; do
  exe=$(basename "$src" .f90)
  gfortran $FFLAGS -I. -o "$exe" $OBJS "$src" -ldl
  (cd "$ROOT" && "$BUILD/$exe")
done
