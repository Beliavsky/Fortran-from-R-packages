#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
BUILD=${1:-"$ROOT/build-gfortran"}
rm -rf "$BUILD"
mkdir -p "$BUILD"
FC=${FC:-gfortran}
FFLAGS=${FFLAGS:-"-std=f2018 -O0 -g -fcheck=all -Wall -Wextra -Werror -Wimplicit-interface"}
COMMON="$FFLAGS -J$BUILD -I$BUILD"

for f in lpsolve_types lpsolve_simplex lpsolve_core lpsolve_special lpsolve; do
  "$FC" $COMMON -c "$ROOT/vendor/lpSolve-fortran-v0.1.0/src/$f.f90" -o "$BUILD/$f.o"
done
for f in adagio_kinds adagio_rng adagio_types adagio_utils adagio_testfunctions \
         adagio_history adagio_geometry adagio_discrete adagio_maxquad \
         adagio_optimize adagio; do
  "$FC" $COMMON -c "$ROOT/src/$f.f90" -o "$BUILD/$f.o"
done
OBJS=$(find "$BUILD" -maxdepth 1 -name '*.o' -printf '%p ')
for src in "$ROOT"/test/*.f90; do
  exe="$BUILD/$(basename "${src%.f90}")"
  "$FC" $COMMON "$src" $OBJS -o "$exe"
  "$exe"
done
