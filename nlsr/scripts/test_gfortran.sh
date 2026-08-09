#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
rm -rf build-strict
mkdir -p build-strict/mod build-strict/obj build-strict/bin
F='-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all'
SRC=(nlsr_kinds nlsr_types nlsr_linalg nlsr_derivatives nlsr_core nlsr_models nlsr_stats nlsr)
for s in "${SRC[@]}"; do
  gfortran $F -Jbuild-strict/mod -Ibuild-strict/mod -c "src/$s.f90" -o "build-strict/obj/$s.o"
done
OBJS=$(printf ' build-strict/obj/%s.o' "${SRC[@]}")
for t in test/*.f90; do
  exe="build-strict/bin/$(basename "$t" .f90)"
  gfortran $F -Ibuild-strict/mod "$t" $OBJS -o "$exe"
  "$exe"
done
for e in example/*.f90; do
  exe="build-strict/bin/$(basename "$e" .f90)"
  gfortran $F -Ibuild-strict/mod "$e" $OBJS -o "$exe"
  "$exe"
done
