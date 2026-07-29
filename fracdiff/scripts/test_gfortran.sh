#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
fc=${FC:-gfortran}

sources=(
  fracdiff_kinds.f90
  fracdiff_status.f90
  fracdiff_types.f90
  fracdiff_linalg.f90
  fracdiff_fft.f90
  fracdiff_difference.f90
  fracdiff_rng.f90
  fracdiff_polynomial.f90
  fracdiff_filter.f90
  fracdiff_optimize.f90
  fracdiff_inference.f90
  fracdiff_simulation.f90
  fracdiff_estimators.f90
  fracdiff_model_api.f90
  fracdiff.f90
)

tests=(
  test_difference
  test_filter
  test_simulation
  test_estimators
  test_fit
)

build_one() {
  local mode=$1
  shift
  local flags=("$@")
  local build="$root/build/$mode"
  local mod="$build/mod"
  local obj="$build/obj"
  rm -rf "$build"
  mkdir -p "$mod" "$obj"

  echo "== $mode library =="
  for source in "${sources[@]}"; do
    "$fc" "${flags[@]}" -J"$mod" -I"$mod" -c "$root/src/$source" \
      -o "$obj/${source%.f90}.o"
  done

  local objects=("$obj"/*.o)

  echo "== $mode tests =="
  "$fc" "${flags[@]}" -J"$mod" -I"$mod" -c "$root/test/test_support.f90" \
    -o "$obj/test_support.o"
  objects+=("$obj/test_support.o")
  for test_name in "${tests[@]}"; do
    "$fc" "${flags[@]}" -J"$mod" -I"$mod" \
      "$root/test/$test_name.f90" "${objects[@]}" -o "$build/$test_name"
    "$build/$test_name"
  done

  echo "== $mode applications and examples =="
  "$fc" "${flags[@]}" -J"$mod" -I"$mod" \
    "$root/app/fracdiff_demo.f90" "$obj"/fracdiff_*.o \
    -o "$build/fracdiff_demo"
  "$build/fracdiff_demo"

  for example_name in long_memory_estimators custom_innovations; do
    "$fc" "${flags[@]}" -J"$mod" -I"$mod" \
      "$root/example/$example_name.f90" "$obj"/fracdiff_*.o \
      -o "$build/$example_name"
    "$build/$example_name"
  done
}

common=(-std=f2018 -Wall -Wextra -Wpedantic -Werror)
build_one strict "${common[@]}" -O0 -g -fcheck=all \
  -fbacktrace -ffpe-trap=invalid,zero,overflow
build_one release "${common[@]}" -O3

echo "All strict and optimized builds passed."
