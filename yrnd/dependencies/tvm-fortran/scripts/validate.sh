#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
build=$(mktemp -d "${TMPDIR:-/tmp}/tvm-fortran.XXXXXX")
trap 'rm -rf "$build"' EXIT

sources=(
  "$root/src/tvm_kinds.f90"
  "$root/src/tvm_root.f90"
  "$root/src/tvm_interpolation.f90"
  "$root/src/tvm_cashflows.f90"
  "$root/src/tvm_curves.f90"
  "$root/src/tvm.f90"
)

flags=(
  -std=f2018
  -Wall
  -Wextra
  -Wpedantic
  -Wconversion-extra
  -Wimplicit-interface
  -Werror
  -fcheck=all
  -fbacktrace
  -O0
)

cd "$build"
gfortran "${flags[@]}" -c "${sources[@]}"
for source in "$root"/test/*.f90; do
  target=$(basename "$source" .f90)
  gfortran "${flags[@]}" ./*.o "$source" -o "$target"
  "./$target"
done
for source in "$root"/app/*.f90 "$root"/example/*.f90; do
  target=$(basename "$source" .f90)
  gfortran "${flags[@]}" ./*.o "$source" -o "$target"
  "./$target" >/dev/null
done

echo "validation: PASS"
