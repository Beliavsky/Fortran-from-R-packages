#!/usr/bin/env sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="$root/build_strict"
flags="-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all"
rm -rf "$build"
mkdir -p "$build/matrix/mod" "$build/piqp/mod" "$build/bin"
msrc="$root/vendor/Matrix-fortran/src"
for s in matrix_kinds matrix_status matrix_dense matrix_decompositions matrix_functions matrix_sparse matrix_sparse_solvers matrix_ordering matrix_io matrix_constructors matrix_sparse_stats matrix_advanced matrix; do
  gfortran $flags -J"$build/matrix/mod" -I"$build/matrix/mod" -c "$msrc/$s.f90" -o "$build/matrix/$s.o"
done
for s in piqp_kinds piqp_types piqp_linalg piqp_solver; do
  gfortran $flags -J"$build/piqp/mod" -I"$build/piqp/mod" -c "$root/src/$s.f90" -o "$build/piqp/$s.o"
done
gfortran $flags -J"$build/piqp/mod" -I"$build/piqp/mod" -I"$build/matrix/mod" -c "$root/src/piqp_matrix_adapter.f90" -o "$build/piqp/piqp_matrix_adapter.o"
gfortran $flags -J"$build/piqp/mod" -I"$build/piqp/mod" -I"$build/matrix/mod" -c "$root/src/piqp.f90" -o "$build/piqp/piqp.o"
for f in "$root"/test/*.f90; do
  n=$(basename "$f" .f90)
  gfortran $flags -I"$build/piqp/mod" -I"$build/matrix/mod" "$f" "$build"/piqp/*.o "$build"/matrix/*.o -o "$build/bin/$n"
  "$build/bin/$n"
done
for f in "$root"/example/*.f90; do
  n=$(basename "$f" .f90)
  gfortran $flags -I"$build/piqp/mod" -I"$build/matrix/mod" "$f" "$build"/piqp/*.o "$build"/matrix/*.o -o "$build/bin/$n"
  "$build/bin/$n"
done
