#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
rm -rf build-strict
mkdir build-strict
FC=${FC:-gfortran}
FLAGS="-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all"
for src in src/trustoptim_kinds.f90 src/trustoptim_types.f90 src/trustoptim_linalg.f90 src/trustoptim.f90 src/trustoptim_binary.f90; do
  "$FC" $FLAGS -J build-strict -I build-strict -c "$src" -o "build-strict/$(basename "${src%.f90}").o"
done
for src in test/*.f90; do
  exe="build-strict/$(basename "${src%.f90}")"
  "$FC" $FLAGS -J build-strict -I build-strict "$src" build-strict/*.o -o "$exe"
  "$exe"
done
for src in example/*.f90; do
  exe="build-strict/$(basename "${src%.f90}")"
  "$FC" $FLAGS -J build-strict -I build-strict "$src" build-strict/*.o -o "$exe"
  "$exe"
done
