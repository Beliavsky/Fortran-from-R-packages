#!/usr/bin/env bash
set -euo pipefail

mode="${1:-debug}"
root="$(cd "$(dirname "$0")" && pwd)"
cd "$root"

fc="${FC:-gfortran}"
common_flags=(-std=f2018 -ffree-line-length-none -Wall -Wextra -Wimplicit-interface -Werror -fbacktrace)
case "$mode" in
  debug)
    mode_flags=(-O0 -g -fcheck=all)
    ;;
  release)
    mode_flags=(-O2)
    ;;
  *)
    echo "usage: ./run_checks.sh [debug|release]" >&2
    exit 2
    ;;
esac

build="build/$mode"
rm -rf "$build"
mkdir -p "$build"

sources=(
  src/msgarch_kinds.f90
  src/msgarch_rng.f90
  src/msgarch_special.f90
  src/msgarch_distributions.f90
  src/msgarch_types.f90
  src/msgarch_models.f90
  src/msgarch_filter.f90
  src/msgarch_simulation.f90
  src/msgarch_optimizer.f90
  src/msgarch_linalg.f90
  src/msgarch_mapping.f90
  src/msgarch_estimation.f90
  src/msgarch_forecast.f90
  src/msgarch_risk.f90
  src/msgarch_hmm.f90
  src/msgarch_posterior.f90
  src/msgarch.f90
)

objects=()
for source in "${sources[@]}"; do
  object="$build/$(basename "${source%.f90}").o"
  "$fc" "${common_flags[@]}" "${mode_flags[@]}" -J"$build" -I"$build" -c "$source" -o "$object"
  objects+=("$object")
done

"$fc" "${common_flags[@]}" "${mode_flags[@]}" -J"$build" -I"$build" -c test/test_helpers.f90 -o "$build/test_helpers.o"

build_program() {
  local source="$1"
  local output="$2"
  "$fc" "${common_flags[@]}" "${mode_flags[@]}" -J"$build" -I"$build" \
    "$source" "$build/test_helpers.o" "${objects[@]}" -o "$build/$output"
}

build_application() {
  local source="$1"
  local output="$2"
  "$fc" "${common_flags[@]}" "${mode_flags[@]}" -J"$build" -I"$build" \
    "$source" "${objects[@]}" -o "$build/$output"
}

build_program test/test_distributions_models.f90 test_distributions_models
build_program test/test_filter_simulation_risk.f90 test_filter_simulation_risk
build_program test/test_estimation_hmm.f90 test_estimation_hmm

build_application app/demo_msgarch.f90 demo_msgarch
build_application example/mcmc_example.f90 mcmc_example
build_application example/fit_csv.f90 fit_csv

"$build/test_distributions_models"
"$build/test_filter_simulation_risk"
"$build/test_estimation_hmm"
"$build/demo_msgarch"
"$build/mcmc_example"
"$build/fit_csv" data/returns.csv single sGARCH norm
"$build/fit_csv" data/returns.csv markov sGARCH norm gjrGARCH std
"$build/fit_csv" data/returns.csv mixture sGARCH norm eGARCH ged

echo "$mode build, tests, and applications passed."
