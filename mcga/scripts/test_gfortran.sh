#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
rm -rf build_strict
mkdir -p build_strict
FC=${FC:-gfortran}
FLAGS=(-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -Jbuild_strict -Ibuild_strict)
SOURCES=(src/mcga_kinds.f90 src/mcga_random.f90 src/mcga_bytes.f90 src/mcga_operators.f90 src/mcga_engine.f90 src/mcga.f90)
OBJS=()
for src in "${SOURCES[@]}"; do
  obj="build_strict/$(basename "${src%.f90}").o"
  "$FC" "${FLAGS[@]}" -c "$src" -o "$obj"
  OBJS+=("$obj")
done
for src in test/*.f90; do
  exe="build_strict/$(basename "${src%.f90}")"
  "$FC" "${FLAGS[@]}" "$src" "${OBJS[@]}" -o "$exe"
  "$exe"
done
for src in example/*.f90; do
  exe="build_strict/$(basename "${src%.f90}")"
  "$FC" "${FLAGS[@]}" "$src" "${OBJS[@]}" -o "$exe"
  "$exe"
done
