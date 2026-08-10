#!/usr/bin/env sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
BUILD="$ROOT/build-strict"
rm -rf "$BUILD"
mkdir -p "$BUILD"
cd "$ROOT"
FLAGS="-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all"
gfortran $FLAGS -J "$BUILD" -I "$BUILD" -c \
  src/ecos_types.f90 src/ecos_sparse.f90 src/ecos_equilibration.f90 src/ecos_linalg.f90 \
  src/ecos_sparse_cones.f90 src/ecos_cones.f90 src/ecos_sparse_solver.f90 \
  src/ecos_solver.f90 src/ecos_bb.f90 src/ecos_api.f90
for src in test/*.f90; do
  exe="$BUILD/$(basename "${src%.f90}")"
  gfortran $FLAGS -J "$BUILD" -I "$BUILD" ./*.o "$src" -o "$exe"
  "$exe"
done
for src in example/*.f90; do
  exe="$BUILD/$(basename "${src%.f90}")"
  gfortran $FLAGS -J "$BUILD" -I "$BUILD" ./*.o "$src" -o "$exe"
  "$exe"
done
rm -f ./*.o ./*.mod
