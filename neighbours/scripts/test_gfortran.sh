#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
BUILD="$ROOT/build_strict"
rm -rf "$BUILD"
mkdir -p "$BUILD"
cd "$BUILD"
FLAGS=(-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all)
gfortran "${FLAGS[@]}" -c \
  "$ROOT/src/neighbours_kinds.f90" \
  "$ROOT/src/neighbours_rng.f90" \
  "$ROOT/src/neighbours.f90"
for src in "$ROOT"/test/*.f90; do
  exe=$(basename "$src" .f90)
  gfortran "${FLAGS[@]}" "$src" neighbours_kinds.o neighbours_rng.o neighbours.o -o "$exe"
  "./$exe"
done
for src in "$ROOT"/example/*.f90; do
  exe=$(basename "$src" .f90)
  gfortran "${FLAGS[@]}" "$src" neighbours_kinds.o neighbours_rng.o neighbours.o -o "$exe"
  "./$exe"
done
