#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"

sources=(
  src/ghyp_kinds.f90
  src/ghyp_special.f90
  src/ghyp_linalg.f90
  src/ghyp_rng.f90
  src/ghyp_gig.f90
  src/ghyp_model.f90
  src/ghyp_distribution.f90
  src/risksimul_types.f90
  src/risksimul_math.f90
  src/risksimul_portfolio.f90
  src/risksimul_simulation.f90
  src/risksimul.f90
)

run_build() {
  local mode="$1"
  shift
  local flags=("$@")
  local dir="build_${mode}"
  rm -rf "$dir" ./*.o ./*.mod
  mkdir -p "$dir"
  gfortran "${flags[@]}" -J "$dir" -I "$dir" -c "${sources[@]}"
  for source in test/*.f90; do
    exe="$dir/$(basename "${source%.f90}")"
    gfortran "${flags[@]}" -J "$dir" -I "$dir" "$source" ./*.o -o "$exe"
    "$exe"
  done
  for source in app/*.f90 example/*.f90; do
    exe="$dir/$(basename "${source%.f90}")"
    gfortran "${flags[@]}" -J "$dir" -I "$dir" "$source" ./*.o -o "$exe"
    "$exe" >/dev/null
  done
}

common=(-std=f2018 -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror)
run_build debug "${common[@]}" -O0 -g -fcheck=all -fbacktrace
run_build optimized "${common[@]}" -O2

echo "validation: PASS"
