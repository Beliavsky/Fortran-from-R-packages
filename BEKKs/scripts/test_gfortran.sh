#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"

fc=${FC:-gfortran}
mode=${1:-all}

sources=(
  src/bekks_kinds.f90
  src/bekks_types.f90
  src/bekks_math.f90
  src/bekks_linalg.f90
  src/bekks_rng.f90
  src/bekks_matrix.f90
  src/bekks_model.f90
  src/bekks_estimation.f90
  src/bekks_forecast.f90
  src/bekks_risk.f90
  src/bekks_diagnostics.f90
  src/bekks_api.f90
)

build_one() {
  local configuration=$1
  local build_dir="build_gfortran/${configuration}"
  local flags
  rm -rf "$build_dir"
  mkdir -p "$build_dir"

  if [[ $configuration == debug ]]; then
    flags=(-std=f2018 -O0 -g -Wall -Wextra -Werror -fcheck=all \
      -ffpe-trap=invalid,zero,overflow -fbacktrace)
  else
    flags=(-std=f2018 -O3 -Wall -Wextra -Werror)
  fi

  local objects=()
  local source object
  for source in "${sources[@]}"; do
    object="$build_dir/$(basename "${source%.f90}").o"
    "$fc" "${flags[@]}" -J"$build_dir" -I"$build_dir" -c "$source" -o "$object"
    objects+=("$object")
  done

  local program executable
  for program in test/*.f90; do
    executable="$build_dir/$(basename "${program%.f90}")"
    "$fc" "${flags[@]}" -J"$build_dir" -I"$build_dir" "$program" \
      "${objects[@]}" -llapack -lblas -o "$executable"
    "$executable"
  done

  for program in app/*.f90 example/*.f90; do
    executable="$build_dir/$(basename "${program%.f90}")"
    "$fc" "${flags[@]}" -J"$build_dir" -I"$build_dir" "$program" \
      "${objects[@]}" -llapack -lblas -o "$executable"
    "$executable" >/dev/null
  done

  echo "$configuration build: PASS"
}

case "$mode" in
  debug|release)
    build_one "$mode"
    ;;
  all)
    build_one debug
    build_one release
    ;;
  *)
    echo "usage: $0 [debug|release|all]" >&2
    exit 2
    ;;
esac
