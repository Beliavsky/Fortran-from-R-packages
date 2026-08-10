#!/usr/bin/env sh
set -eu
cd "$(dirname "$0")/.."
rm -rf build-strict
mkdir build-strict
flags="-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all"
gfortran $flags -J build-strict -I build-strict -c src/nnls_kinds.f90 -o build-strict/nnls_kinds.o
gfortran $flags -J build-strict -I build-strict -c src/nnls_linalg.f90 -o build-strict/nnls_linalg.o
gfortran $flags -J build-strict -I build-strict -c src/nnls_solver.f90 -o build-strict/nnls_solver.o
gfortran $flags -J build-strict -I build-strict -c src/nnls.f90 -o build-strict/nnls.o
objs="build-strict/nnls_kinds.o build-strict/nnls_linalg.o build-strict/nnls_solver.o build-strict/nnls.o"
for src in test/*.f90; do
  exe="build-strict/$(basename "${src%.f90}")"
  gfortran $flags -I build-strict $objs "$src" -o "$exe"
  "$exe"
done
for src in example/*.f90; do
  exe="build-strict/$(basename "${src%.f90}")"
  gfortran $flags -I build-strict $objs "$src" -o "$exe"
  "$exe"
done
