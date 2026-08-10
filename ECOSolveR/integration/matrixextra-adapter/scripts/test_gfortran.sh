#!/usr/bin/env sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
ECOS=$(CDPATH= cd -- "$ROOT/../.." && pwd)
MX="$ROOT/vendor/MatrixExtra-fortran"
BUILD="$ROOT/build-strict"
FLAGS="-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all"
rm -rf "$BUILD"
mkdir -p "$BUILD/mod" "$BUILD/obj" "$BUILD/bin"

for s in matrix_kinds matrix_status matrix_dense matrix_decompositions matrix_functions matrix_sparse \
         matrix_sparse_solvers matrix_ordering matrix_io matrix_constructors matrix_sparse_stats \
         matrix_advanced matrix; do
  gfortran $FLAGS -J "$BUILD/mod" -I "$BUILD/mod" -c \
    "$MX/vendor/Matrix-fortran/src/$s.f90" -o "$BUILD/obj/$s.o"
done
for s in matrixextra_types matrixextra_conversions matrixextra_utils matrixextra_slice matrixextra_bind \
         matrixextra_matmul matrixextra_ops matrixextra_linalg matrixextra_recycle matrixextra_pattern matrixextra; do
  gfortran $FLAGS -J "$BUILD/mod" -I "$BUILD/mod" -c "$MX/src/$s.f90" -o "$BUILD/obj/$s.o"
done
for s in ecos_types ecos_sparse ecos_equilibration ecos_linalg ecos_sparse_cones ecos_cones ecos_sparse_solver ecos_solver ecos_bb ecos_api; do
  gfortran $FLAGS -J "$BUILD/mod" -I "$BUILD/mod" -c "$ECOS/src/$s.f90" -o "$BUILD/obj/$s.o"
done
gfortran $FLAGS -J "$BUILD/mod" -I "$BUILD/mod" -c "$ROOT/src/ecos_matrixextra_adapter.f90" \
  -o "$BUILD/obj/ecos_matrixextra_adapter.o"
gfortran $FLAGS -I "$BUILD/mod" "$ROOT/example/matrixextra_sparse_lp.f90" "$BUILD"/obj/*.o \
  -o "$BUILD/bin/matrixextra_sparse_lp"
"$BUILD/bin/matrixextra_sparse_lp"
