#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"

mode=${1:-checked}
case "$mode" in
  checked)
    flags=(-std=f2018 -O0 -Wall -Wextra -Wconversion-extra -Werror \
      -fcheck=all -fbacktrace)
    ;;
  optimized)
    flags=(-std=f2018 -O2 -Wall -Wextra -Wconversion-extra -Werror)
    ;;
  *)
    echo "usage: $0 [checked|optimized]" >&2
    exit 2
    ;;
esac

sources=(
  src/ghyp_kinds.f90
  src/ghyp_special.f90
  src/ghyp_rng.f90
  src/ghyp_linalg.f90
  src/ghyp_gig.f90
  src/ghyp_model.f90
  src/ghyp_distribution.f90
  src/ghyp_risk.f90
  src/ghyp_optimize.f90
  src/ghyp_fitting.f90
  src/ghyp_portfolio.f90
  src/ghyp_utilities.f90
  src/ghyp.f90
  src/sharpe_rratio_calibration.f90
  src/sharpe_rratio_records.f90
  src/sharpe_rratio_statistics.f90
  src/sharpe_rratio_estimator.f90
  src/sharpe_rratio.f90
)

build="build-$mode"
rm -rf "$build"
mkdir -p "$build/mod" "$build/obj" "$build/bin"

for source in "${sources[@]}"; do
  object="$build/obj/$(basename "${source%.f90}").o"
  gfortran "${flags[@]}" -J "$build/mod" -I "$build/mod" \
    -c "$source" -o "$object"
done

for source in test/*.f90 app/*.f90 example/*.f90; do
  name=$(basename "$source" .f90)
  gfortran "${flags[@]}" -I "$build/mod" "$source" "$build"/obj/*.o \
    -o "$build/bin/$name"
done
