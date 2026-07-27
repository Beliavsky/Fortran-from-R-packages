#!/usr/bin/env bash
set -euo pipefail

fc=${FC:-gfortran}
root=$(cd "$(dirname "$0")/.." && pwd)

sources=(
  greeks_kinds.f90
  greeks_types.f90
  greeks_math.f90
  greeks_rng.f90
  greeks_european.f90
  greeks_geometric_asian.f90
  greeks_american.f90
  greeks_paths.f90
  greeks_malliavin.f90
  greeks_implied_volatility.f90
  greeks.f90
)

run_build() {
  local mode=$1
  shift
  local flags=("$@")
  local build="$root/build/$mode"
  rm -rf "$build"
  mkdir -p "$build"
  cd "$build"
  local source_paths=()
  local source
  for source in "${sources[@]}"; do
    source_paths+=("$root/src/$source")
  done
  "$fc" "${flags[@]}" -J. -I. -c "${source_paths[@]}"
  local file name
  for file in "$root"/test/*.f90; do
    name=$(basename "$file" .f90)
    "$fc" "${flags[@]}" -J. -I. ./*.o "$file" -o "$name"
    "./$name"
  done
  for file in "$root"/app/*.f90 "$root"/example/*.f90; do
    name=$(basename "$file" .f90)
    "$fc" "${flags[@]}" -J. -I. ./*.o "$file" -o "$name"
    "./$name" >/dev/null
  done
}

common=(-std=f2018 -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror)
run_build debug "${common[@]}" -O0 -fcheck=all -fbacktrace
run_build release "${common[@]}" -O2

echo "validation: PASS"
