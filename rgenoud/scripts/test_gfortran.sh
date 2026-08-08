#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
flags=(-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all)
rm -rf build_gfortran
mkdir build_gfortran
sources=(
  src/rgenoud_kinds.f90
  src/rgenoud_random.f90
  src/rgenoud_types.f90
  src/rgenoud_derivatives.f90
  src/rgenoud_operators.f90
  src/rgenoud_core.f90
  src/rgenoud_stats.f90
  src/rgenoud.f90
)
objects=()
for src in "${sources[@]}"; do
  obj="build_gfortran/$(basename "${src%.f90}").o"
  gfortran "${flags[@]}" -J build_gfortran -I build_gfortran -c "$src" -o "$obj"
  objects+=("$obj")
done
for test_src in test/*.f90; do
  exe="build_gfortran/$(basename "${test_src%.f90}")"
  gfortran "${flags[@]}" -J build_gfortran -I build_gfortran \
    "$test_src" "${objects[@]}" -o "$exe"
  "$exe"
done
