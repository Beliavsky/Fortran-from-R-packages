#!/usr/bin/env sh
set -eu
cd "$(dirname "$0")/.."
flags="-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all"
build=build_strict
rm -rf "$build"
mkdir -p "$build/mod" "$build/obj" "$build/bin"
vendor="vendor/Matrix-fortran"
for s in matrix_kinds matrix_status matrix_dense matrix_decompositions matrix_functions matrix_sparse matrix_sparse_solvers matrix_ordering matrix_io matrix_constructors matrix_sparse_stats matrix_advanced matrix; do
  gfortran $flags -J "$build/mod" -I "$build/mod" -c "$vendor/src/$s.f90" -o "$build/obj/$s.o"
done
for s in matrixextra_types matrixextra_conversions matrixextra_utils matrixextra_slice matrixextra_bind matrixextra_matmul matrixextra_ops matrixextra_linalg matrixextra_recycle matrixextra_pattern matrixextra; do
  gfortran $flags -J "$build/mod" -I "$build/mod" -c "src/$s.f90" -o "$build/obj/$s.o"
done
for src in test/*.f90 example/*.f90; do
  name=$(basename "$src" .f90)
  gfortran $flags -I "$build/mod" "$src" "$build"/obj/*.o -o "$build/bin/$name"
done
for p in "$build"/bin/test_*; do "$p"; done
for p in "$build"/bin/sparse_pipeline "$build"/bin/matmul_example; do "$p"; done
