#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

modules=(
  qrmtools_kinds qrmtools_types qrmtools_stats qrmtools_optimization
  qrmtools_distributions qrmtools_evt qrmtools_brownian
  qrmtools_black_scholes qrmtools_risk qrmtools_bounds qrmtools_returns
  qrmtools_hierarchy qrmtools_allocation qrmtools_tests qrmtools_garch
  qrmtools
)

build_one() {
  local build_dir="$1"
  local optimization="$2"
  rm -rf "$build_dir"
  mkdir -p "$build_dir"

  local flags=(
    -std=f2018 "$optimization" -Wall -Wextra -Wconversion-extra
    -Wimplicit-interface -Werror -fcheck=all -fbacktrace
    -J "$build_dir" -I "$build_dir"
  )

  local module
  for module in "${modules[@]}"; do
    gfortran "${flags[@]}" -c "src/${module}.f90" -o "$build_dir/${module}.o"
  done

  local objects=()
  for module in "${modules[@]}"; do
    objects+=("$build_dir/${module}.o")
  done

  local source executable
  for source in test/*.f90 app/*.f90 example/*.f90; do
    executable="$build_dir/$(basename "${source%.f90}")"
    gfortran "${flags[@]}" "$source" "${objects[@]}" -o "$executable"
    "$executable"
  done
}

build_one build-validation-debug -O0
build_one build-validation-optimized -O2

echo "validation: PASS"
