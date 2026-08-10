#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
rm -rf build_strict
mkdir build_strict
FLAGS=(-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all)
SRC=(cgnm_kinds cgnm_types cgnm_utils cgnm_linalg cgnm_kmeans cgnm_core cgnm_postprocess cgnm_extensions cgnm)
for s in "${SRC[@]}"; do
  gfortran "${FLAGS[@]}" -J build_strict -I build_strict -c "src/$s.f90" -o "build_strict/$s.o"
done
OBJ=(build_strict/*.o)
for t in test/*.f90; do
  exe="build_strict/$(basename "${t%.f90}")"
  gfortran "${FLAGS[@]}" -I build_strict "$t" "${OBJ[@]}" -o "$exe"
  "$exe"
done
for e in example/*.f90; do
  exe="build_strict/$(basename "${e%.f90}")"
  gfortran "${FLAGS[@]}" -I build_strict "$e" "${OBJ[@]}" -o "$exe"
  "$exe"
done
