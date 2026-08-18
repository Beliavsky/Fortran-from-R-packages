#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
build="$root/build-optimized"
rm -rf "$build"
mkdir -p "$build"
cd "$build"
flags=(-std=f2018 -O2 -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror)
sources=(
  "$root/src/copula_kinds.f90"
  "$root/src/copula_types.f90"
  "$root/src/copula_special.f90"
  "$root/src/copula_random.f90"
  "$root/src/copula_linalg.f90"
  "$root/src/copula_families.f90"
  "$root/src/copula_simulation.f90"
  "$root/src/copula_dependence.f90"
  "$root/src/copula_empirical.f90"
  "$root/src/copula_fitting.f90"
  "$root/src/copula_compositions.f90"
  "$root/src/copula_special_discrete.f90"
  "$root/src/copula.f90"
)
gfortran "${flags[@]}" -J . -I . -c "${sources[@]}"
for source in "$root"/test/*.f90; do
  exe=$(basename "$source" .f90)
  gfortran "${flags[@]}" -J . -I . "$source" ./*.o -o "$exe"
  "./$exe"
done
echo 'optimized validation: PASS'
