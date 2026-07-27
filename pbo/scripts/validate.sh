#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
rm -rf build-validation
mkdir build-validation
modflags=(-Jbuild-validation -Ibuild-validation)
debug=(-std=f2018 -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror -fcheck=all -fbacktrace)
sources=(
  src/pbo_kinds.f90
  src/pbo_types.f90
  src/pbo_combinations.f90
  src/pbo_stats.f90
  src/pbo_metrics.f90
  src/pbo_core.f90
  src/pbo_analysis.f90
  src/pbo.f90
)
objects=()
for source in "${sources[@]}"; do
  object="build-validation/$(basename "${source%.f90}").o"
  gfortran "${debug[@]}" "${modflags[@]}" -c "$source" -o "$object"
  objects+=("$object")
done
for source in test/*.f90 app/*.f90 example/*.f90; do
  exe="build-validation/$(basename "${source%.f90}")"
  gfortran "${debug[@]}" "${modflags[@]}" "$source" "${objects[@]}" -o "$exe"
  "$exe"
done
printf '%s\n' 'validation: PASS'
