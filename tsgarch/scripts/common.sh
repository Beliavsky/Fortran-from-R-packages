#!/usr/bin/env sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
FC=${FC:-gfortran}
BUILD=${BUILD:-$ROOT/build}
MOD="$BUILD/mod"
OBJ="$BUILD/obj"
BIN="$BUILD/bin"
mkdir -p "$MOD" "$OBJ" "$BIN"
SOURCES="
$ROOT/src/ghyp_kinds.f90
$ROOT/src/ghyp_special.f90
$ROOT/src/ghyp_rng.f90
$ROOT/src/ghyp_linalg.f90
$ROOT/src/ghyp_gig.f90
$ROOT/src/ghyp_model.f90
$ROOT/src/ghyp_distribution.f90
$ROOT/src/tsd_types.f90
$ROOT/src/tsd_math.f90
$ROOT/src/tsd_optimize.f90
$ROOT/src/tsd_distributions.f90
$ROOT/src/tsd_fit.f90
$ROOT/src/tsd_moments.f90
$ROOT/src/tsd_spd.f90
$ROOT/src/tsd_profile.f90
$ROOT/src/tsdistributions.f90
$ROOT/src/tsgarch_types.f90
$ROOT/src/tsgarch_model.f90
$ROOT/src/tsgarch_fit.f90
$ROOT/src/tsgarch_simulation.f90
$ROOT/src/tsgarch_forecast.f90
$ROOT/src/tsgarch_diagnostics.f90
$ROOT/src/tsgarch_profile.f90
$ROOT/src/tsgarch_backtest.f90
$ROOT/src/tsgarch_benchmarks.f90
$ROOT/src/tsgarch.f90
"
TESTS="test_models test_simulation_forecast test_estimation test_inference_profile test_backtest test_constraints_vreg"
compile_library() {
  rm -f "$OBJ"/*.o "$MOD"/*.mod
  # shellcheck disable=SC2086
  "$FC" $FLAGS -J"$MOD" -I"$MOD" -c $SOURCES
  mv ./*.o "$OBJ"/
  "$FC" $FLAGS -J"$MOD" -I"$MOD" -c "$ROOT/test/test_support.f90" -o "$OBJ/test_support.o"
}
run_tests() {
  for name in $TESTS; do
    "$FC" $FLAGS -J"$MOD" -I"$MOD" "$ROOT/test/$name.f90" "$OBJ"/*.o -o "$BIN/$name"
    "$BIN/$name"
  done
}
