#!/usr/bin/env bash
set -euo pipefail
mode="${1:-check}"
root="$(cd "$(dirname "$0")/.." && pwd)"
build="$root/build/$mode"
mod="$build/mod"
mkdir -p "$build" "$mod"
common=(-std=f2018 -ffree-line-length-none -Wall -Wextra -Werror -Wno-maybe-uninitialized -J"$mod" -I"$mod")
if [[ "$mode" == "check" ]]; then
  flags=("${common[@]}" -O0 -g -fcheck=all -fbacktrace -finit-real=snan -finit-integer=-999999)
else
  flags=("${common[@]}" -O3 -DNDEBUG)
fi
sources=(
  tscopula_kinds tscopula_status tscopula_math tscopula_margins
  tscopula_vtransforms tscopula_paircopula tscopula_timeseries
  tscopula_dvine tscopula_models tscopula_compat tscopula
)
objects=()
for name in "${sources[@]}"; do
  obj="$build/$name.o"
  gfortran "${flags[@]}" -c "$root/src/$name.f90" -o "$obj"
  objects+=("$obj")
done
gfortran "${flags[@]}" -c "$root/test/test_utils.f90" -o "$build/test_utils.o"
objects+=("$build/test_utils.o")
tests=(test_margins_vtrans test_paircopula test_timeseries test_dvine test_models test_fitting)
for name in "${tests[@]}"; do
  exe="$build/$name"
  gfortran "${flags[@]}" "$root/test/$name.f90" "${objects[@]}" -o "$exe"
  "$exe"
done
