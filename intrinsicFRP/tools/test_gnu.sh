#!/usr/bin/env bash
set -euo pipefail

mode="${1:-debug}"
root="$(cd "$(dirname "$0")/.." && pwd)"
build="$root/build/$mode"
rm -rf "$build"
mkdir -p "$build/mod" "$build/obj" "$build/bin"

common=(-std=f2018 -Wall -Wextra -Wpedantic -Werror \
  -Wno-maybe-uninitialized -ffree-line-length-132)
if [[ "$mode" == "debug" ]]; then
  flags=("${common[@]}" -O0 -g -fcheck=all -fbacktrace \
    -ffpe-trap=invalid,zero,overflow)
elif [[ "$mode" == "release" ]]; then
  flags=("${common[@]}" -O3)
else
  echo "usage: $0 [debug|release]" >&2
  exit 2
fi

sources=(
  intrinsicfrp_kinds.f90
  intrinsicfrp_types.f90
  intrinsicfrp_linalg.f90
  intrinsicfrp_stats.f90
  intrinsicfrp_hac.f90
  intrinsicfrp_models.f90
  intrinsicfrp_identification.f90
  intrinsicfrp_oracle.f90
  intrinsicfrp.f90
)

for source in "${sources[@]}"; do
  object="$build/obj/${source%.f90}.o"
  gfortran "${flags[@]}" -J"$build/mod" -I"$build/mod" \
    -c "$root/src/$source" -o "$object"
done
ar rcs "$build/libintrinsicfrp.a" "$build"/obj/*.o

run_group() {
  local directory="$1"
  local source name executable
  shopt -s nullglob
  for source in "$root/$directory"/*.f90; do
    name="$(basename "${source%.f90}")"
    executable="$build/bin/$name"
    gfortran "${flags[@]}" -I"$build/mod" "$source" \
      "$build/libintrinsicfrp.a" -o "$executable"
    "$executable"
  done
}

run_group test
run_group example
run_group app
printf 'GNU Fortran %s validation: PASS\n' "$mode"
