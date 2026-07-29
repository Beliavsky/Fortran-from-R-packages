#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"

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
)

run_build() {
  local mode=$1
  shift
  local dir="build-validation-$mode"
  rm -rf "$dir"
  mkdir -p "$dir/mod" "$dir/obj" "$dir/bin"
  gfortran "$@" -J "$dir/mod" -I "$dir/mod" -c "${sources[@]}"
  mv ./*.o "$dir/obj/"
  for source in test/*.f90 app/*.f90 example/*.f90; do
    name=$(basename "$source" .f90)
    gfortran "$@" -I "$dir/mod" "$source" "$dir"/obj/*.o -o "$dir/bin/$name"
  done
  for executable in "$dir"/bin/test_*; do
    "$executable"
  done
  "$dir/bin/ghyp_demo" >/dev/null
  "$dir/bin/distributions_and_risk" >/dev/null
  "$dir/bin/fitting_and_portfolio" >/dev/null
}

run_build debug -std=f2018 -O0 -Wall -Wextra -Wconversion-extra -Werror -fcheck=all -fbacktrace
run_build optimized -std=f2018 -O2 -Wall -Wextra -Wconversion-extra -Werror

echo "validation: PASS"
