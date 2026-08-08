#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
rm -rf build-gfortran
mkdir -p build-gfortran/mod build-gfortran/obj build-gfortran/bin
flags=(-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -Jbuild-gfortran/mod -Ibuild-gfortran/mod)
sources=(
  src/pso_kinds.f90
  src/pso_types.f90
  src/pso_random.f90
  src/pso_lbfgsb.f90
  src/pso_core.f90
  src/pso_benchmarks.f90
  src/pso.f90
)
objects=()
for src in "${sources[@]}"; do
  obj="build-gfortran/obj/$(basename "${src%.f90}").o"
  gfortran "${flags[@]}" -c "$src" -o "$obj"
  objects+=("$obj")
done
for src in test/*.f90; do
  exe="build-gfortran/bin/$(basename "${src%.f90}")"
  gfortran "${flags[@]}" "$src" "${objects[@]}" -o "$exe"
  "$exe"
done
for src in example/*.f90; do
  exe="build-gfortran/bin/$(basename "${src%.f90}")"
  gfortran "${flags[@]}" "$src" "${objects[@]}" -o "$exe"
  "$exe"
done
