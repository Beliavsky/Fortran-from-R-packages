#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
build="$root/build-gfortran"
rm -rf "$build" && mkdir -p "$build"
cd "$build"
flags=(-std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow)
src=(imputefin_kinds imputefin_types imputefin_rng imputefin_math imputefin_linalg imputefin_missing imputefin_ar1_gaussian imputefin_ar1_t imputefin_var_t imputefin_wrappers imputefin)
for s in "${src[@]}"; do gfortran "${flags[@]}" -J. -I. -c "$root/src/$s.f90"; done
objs=(); for s in "${src[@]}"; do objs+=("$s.o"); done
for t in "$root"/test/*.f90; do
  exe=$(basename "$t" .f90)
  gfortran "${flags[@]}" -I. "$t" "${objs[@]}" -o "$exe"
  "./$exe"
done
