#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
build="$root/build-gfortran"
rm -rf "$build"
mkdir -p "$build/mod" "$build/obj" "$build/bin"
flags=(-std=f2018 -Wall -Wextra -Wpedantic -Wimplicit-interface -fimplicit-none -J"$build/mod" -I"$build/mod")
if [[ "${DEBUG:-1}" == "1" ]]; then
  flags+=(-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace)
else
  flags+=(-O3)
fi
sources=(
  yieldcurve_kinds.f90
  yieldcurve_status.f90
  yieldcurve_factors.f90
  yieldcurve_linalg.f90
  yieldcurve_optimization.f90
  yieldcurve_models.f90
  yieldcurve.f90
)
objects=()
for source in "${sources[@]}"; do
  object="$build/obj/${source%.f90}.o"
  gfortran "${flags[@]}" -c "$root/src/$source" -o "$object"
  objects+=("$object")
done
gfortran "${flags[@]}" "$root/test/test_yieldcurve.f90" "${objects[@]}" -o "$build/bin/test_yieldcurve"
gfortran "${flags[@]}" "$root/example/basic_fit.f90" "${objects[@]}" -o "$build/bin/basic_fit"
"$build/bin/test_yieldcurve"
"$build/bin/basic_fit"
