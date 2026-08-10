#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
build="$root/build-strict"
rm -rf "$build"
mkdir -p "$build/mod"
cd "$build"
flags=(-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all)
sources=(bvls_kinds.f90 bvls_types.f90 bvls_qr.f90 bvls_solver.f90 bvls.f90)
objects=()
for src in "${sources[@]}"; do
  obj="${src%.f90}.o"
  gfortran "${flags[@]}" -Jmod -Imod -c "$root/src/$src" -o "$obj"
  objects+=("$obj")
done
for t in "$root"/test/*.f90; do
  exe="$(basename "${t%.f90}")"
  gfortran "${flags[@]}" -Imod "$t" "${objects[@]}" -o "$exe"
  "./$exe"
done
for e in "$root"/example/*.f90; do
  exe="$(basename "${e%.f90}")"
  gfortran "${flags[@]}" -Imod "$e" "${objects[@]}" -o "$exe"
  "./$exe"
done
