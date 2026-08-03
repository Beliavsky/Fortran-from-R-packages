#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
rm -rf build-gfortran
mkdir -p build-gfortran/mod build-gfortran/obj build-gfortran/bin
flags=(-std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all \
  -ffpe-trap=invalid,zero,overflow -O0 -g)
sources=(
  src/nloptr_kinds.f90
  src/nloptr_types.f90
  src/nloptr_utils.f90
  src/nloptr_derivatives.f90
  src/nloptr_evaluation.f90
  src/nloptr_solvers.f90
  src/nloptr_api.f90
  src/nloptr.f90
  src/nloptr_example_functions.f90
)
for source in "${sources[@]}"; do
  object="build-gfortran/obj/$(basename "${source%.f90}").o"
  gfortran "${flags[@]}" -J build-gfortran/mod -I build-gfortran/mod -c "$source" -o "$object"
done
objects=(build-gfortran/obj/*.o)
for source in test/*.f90 example/*.f90 app/*.f90; do
  executable="build-gfortran/bin/$(basename "${source%.f90}")"
  gfortran "${flags[@]}" -J build-gfortran/mod -I build-gfortran/mod \
    "$source" "${objects[@]}" -o "$executable"
  "$executable"
done
