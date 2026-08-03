#!/usr/bin/env sh
set -eu
cd "$(dirname "$0")/.."
SOURCES="src/optimx_kinds.f90 src/optimx_types.f90 src/optimx_linalg.f90 \
src/optimx_eval.f90 src/optimx_solvers.f90 src/optimx_checks.f90 \
src/optimx_api.f90 src/optimx.f90 src/optimx_example_functions.f90"
run_build() {
  name="$1"; shift
  rm -rf "build-$name"
  mkdir "build-$name"
  gfortran "$@" -J "build-$name" -I "build-$name" -c $SOURCES
  for file in test/*.f90; do
    exe="build-$name/$(basename "${file%.f90}")"
    gfortran "$@" -J "build-$name" -I "build-$name" "$file" ./*.o -o "$exe"
    "$exe"
  done
  for file in example/*.f90 app/*.f90; do
    exe="build-$name/$(basename "${file%.f90}")"
    gfortran "$@" -J "build-$name" -I "build-$name" "$file" ./*.o -o "$exe"
    "$exe" >/dev/null
  done
  rm -f ./*.o
}
run_build debug -std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all \
  -ffpe-trap=invalid,zero,overflow -fbacktrace
run_build release -std=f2018 -O3 -Wall -Wextra -Wpedantic -Werror
