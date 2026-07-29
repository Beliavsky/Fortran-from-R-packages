#!/usr/bin/env bash
set -euo pipefail

mode="${1:-strict}"
fc="${FC:-gfortran}"
root="$(cd "$(dirname "$0")" && pwd)"
build="$root/build/$mode"
rm -rf "$build"
mkdir -p "$build/mod" "$build/obj" "$build/bin"

common=(-std=f2018 -ffree-line-length-none -Wall -Wextra -Werror -Wimplicit-interface)
case "$mode" in
  strict)
    flags=("${common[@]}" -Wconversion-extra -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow)
    ;;
  optimized)
    flags=("${common[@]}" -O3)
    ;;
  *)
    echo "usage: $0 [strict|optimized]" >&2
    exit 2
    ;;
esac

sources=(
  acdm_kinds.f90
  acdm_math.f90
  acdm_distributions.f90
  acdm_models.f90
  acdm_fit.f90
  acdm_data.f90
  acdm_diagnostics.f90
  acdm_profiles.f90
  acdm.f90
  acdm_api.f90
)
objects=()
for src in "${sources[@]}"; do
  obj="$build/obj/${src%.f90}.o"
  "$fc" "${flags[@]}" -J "$build/mod" -I "$build/mod" \
    -c "$root/src/$src" -o "$obj"
  objects+=("$obj")
done

build_and_run_dir() {
  local directory="$1"
  local file name exe
  shopt -s nullglob
  for file in "$root/$directory"/*.f90; do
    name="$(basename "${file%.f90}")"
    exe="$build/bin/$name"
    "$fc" "${flags[@]}" -I "$build/mod" "$file" "${objects[@]}" -o "$exe"
    echo "==> $directory/$name"
    "$exe"
  done
  shopt -u nullglob
}

build_and_run_dir test
build_and_run_dir app
build_and_run_dir example

echo "All $mode targets passed."
