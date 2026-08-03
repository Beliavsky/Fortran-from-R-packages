#!/usr/bin/env bash
# SPDX-License-Identifier: Artistic-2.0
set -euo pipefail

mode="${1:-strict}"
root="$(cd "$(dirname "$0")/.." && pwd)"
build="$root/build/gfortran-$mode"
mods="$build/mod"
objs="$build/obj"
bin="$build/bin"

case "$mode" in
  strict)
    flags=(-std=f2018 -O0 -g -Wall -Wextra -Werror -Wno-compare-reals \
      -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace)
    ;;
  release)
    flags=(-std=f2018 -O3 -Wall -Wextra -Werror -Wno-compare-reals)
    ;;
  *)
    echo "usage: $0 [strict|release]" >&2
    exit 2
    ;;
esac

rm -rf "$build"
mkdir -p "$mods" "$objs" "$bin"

sources=(
  ecd_kinds
  ecd_rng
  ecd_math
  ecd_core
  ecld_models
  ecd_processes
  ecd_options
  ecd_timeseries
  lamp_process
  ecd_fitting
  ecd_compat
  ecd_api
)

object_files=()
for name in "${sources[@]}"; do
  object="$objs/$name.o"
  gfortran "${flags[@]}" -J"$mods" -I"$mods" -c \
    "$root/src/$name.f90" -o "$object"
  object_files+=("$object")
done

run_group() {
  local group="$1"
  local source exe name
  shopt -s nullglob
  for source in "$root/$group"/*.f90; do
    name="$(basename "${source%.f90}")"
    exe="$bin/$name"
    gfortran "${flags[@]}" -I"$mods" "$source" \
      "${object_files[@]}" -o "$exe"
    echo "[$mode] $group/$name"
    "$exe"
  done
}

run_group test
run_group app
run_group example

echo "All $mode builds and runs passed."
