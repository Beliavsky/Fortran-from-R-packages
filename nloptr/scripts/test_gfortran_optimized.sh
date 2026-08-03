#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
rm -rf build-gfortran-opt
mkdir -p build-gfortran-opt/mod build-gfortran-opt/obj build-gfortran-opt/bin
flags=(-std=f2018 -Wall -Wextra -Wpedantic -Werror -O3)
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
  object="build-gfortran-opt/obj/$(basename "${source%.f90}").o"
  gfortran "${flags[@]}" -J build-gfortran-opt/mod -I build-gfortran-opt/mod -c "$source" -o "$object"
done
objects=(build-gfortran-opt/obj/*.o)
for source in test/*.f90 example/*.f90 app/*.f90; do
  executable="build-gfortran-opt/bin/$(basename "${source%.f90}")"
  gfortran "${flags[@]}" -J build-gfortran-opt/mod -I build-gfortran-opt/mod \
    "$source" "${objects[@]}" -o "$executable"
  "$executable"
done
