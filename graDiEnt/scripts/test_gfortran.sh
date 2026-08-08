#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="$ROOT/build-strict"
rm -rf "$BUILD"
mkdir -p "$BUILD/mod" "$BUILD/obj" "$BUILD/bin"
FLAGS=(-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -J"$BUILD/mod" -I"$BUILD/mod")
SRC=(gradient_kinds gradient_rng gradient_types gradient_stats gradient_sqgde gradient_benchmarks gradient)
for s in "${SRC[@]}"; do
  gfortran "${FLAGS[@]}" -c "$ROOT/src/$s.f90" -o "$BUILD/obj/$s.o"
done
OBJS=("$BUILD"/obj/*.o)
for t in "$ROOT"/test/*.f90; do
  name="$(basename "$t" .f90)"
  gfortran "${FLAGS[@]}" "$t" "${OBJS[@]}" -o "$BUILD/bin/$name"
  "$BUILD/bin/$name"
done
for e in "$ROOT"/example/*.f90; do
  name="$(basename "$e" .f90)"
  gfortran "${FLAGS[@]}" "$e" "${OBJS[@]}" -o "$BUILD/bin/$name"
  "$BUILD/bin/$name"
done
