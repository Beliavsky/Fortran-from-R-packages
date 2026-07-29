#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
fc=${FC:-gfortran}

sources=(
  kind_mod.f90
  stats_mod.f90
  special_functions_mod.f90
  linalg_mod.f90
  optimization_mod.f90
  random_mod.f90
  distribution_mod.f90
  rugarch.f90
  rugarch_diagnostics.f90
  rugarch_extensions.f90
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
  smoots_kinds.f90
  smoots_status.f90
  smoots_types.f90
  smoots_stats.f90
  smoots_linalg.f90
  smoots_smoothing.f90
  smoots_arma.f90
  smoots_estimation.f90
  smoots_forecast.f90
  smoots.f90
  ufrisk_types.f90
  ufrisk_math.f90
  ufrisk_backtests.f90
  ufrisk_smoothing.f90
  ufrisk_loggarch.f90
  ufrisk_varcast.f90
  ufrisk.f90
)

tests=(
  test_core_algorithms
  test_traffic
  test_varcast_models
  test_semiparametric
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

  local library_objects=("$obj"/*.o)
  "$fc" "${flags[@]}" -J"$mod" -I"$mod" -c "$root/test/test_support.f90" \
    -o "$obj/test_support.o"

  echo "== $mode tests =="
  for test_name in "${tests[@]}"; do
    "$fc" "${flags[@]}" -J"$mod" -I"$mod" "$root/test/$test_name.f90" \
      "${library_objects[@]}" "$obj/test_support.o" -o "$build/$test_name"
    "$build/$test_name"
  done

  echo "== $mode applications and examples =="
  "$fc" "${flags[@]}" -J"$mod" -I"$mod" "$root/app/ufrisk_demo.f90" \
    "${library_objects[@]}" -o "$build/ufrisk_demo"
  "$build/ufrisk_demo"
  "$fc" "${flags[@]}" -J"$mod" -I"$mod" "$root/example/semiparametric_figarch.f90" \
    "${library_objects[@]}" -o "$build/semiparametric_figarch"
  "$build/semiparametric_figarch"
}

common=(-std=f2018 -Wall -Wextra -Wpedantic -Werror -fimplicit-none -ffree-line-length-none)
build_one strict "${common[@]}" -O0 -g -fcheck=all -fbacktrace \
  -ffpe-trap=invalid,zero,overflow
build_one release "${common[@]}" -O3

echo "All strict and optimized builds passed."
