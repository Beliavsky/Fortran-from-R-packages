#!/usr/bin/env sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
BUILD=${TMPDIR:-/tmp}/irlba-fortran-strict
rm -rf "$BUILD"
mkdir -p "$BUILD"
FLAGS="-std=f2018 -O2 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all"
for f in irlba_kinds irlba_lapack irlba_sparse irlba_linalg irlba_operator irlba_core irlba_algorithms irlba; do
  gfortran $FLAGS -J"$BUILD" -I"$BUILD" -c "$ROOT/src/$f.f90" -o "$BUILD/$f.o"
done
OBJS="$BUILD/irlba_kinds.o $BUILD/irlba_lapack.o $BUILD/irlba_sparse.o $BUILD/irlba_linalg.o $BUILD/irlba_operator.o $BUILD/irlba_core.o $BUILD/irlba_algorithms.o $BUILD/irlba.o"
for src in "$ROOT"/test/*.f90 "$ROOT"/example/basic.f90; do
  exe="$BUILD/$(basename "${src%.f90}")"
  gfortran $FLAGS -J"$BUILD" -I"$BUILD" "$src" $OBJS -llapack -lblas -o "$exe"
  "$exe"
done
