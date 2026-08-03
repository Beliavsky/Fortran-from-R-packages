#!/usr/bin/env bash
set -euo pipefail
FC=${FC:-gfortran}
ROOT=$(cd "$(dirname "$0")" && pwd)
BUILD="$ROOT/build-direct"
rm -rf "$BUILD"
mkdir -p "$BUILD/mod" "$BUILD/obj" "$BUILD/bin"
FLAGS=(-std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all \
  -ffpe-trap=invalid,zero,overflow -fbacktrace -O0 -g \
  -J"$BUILD/mod" -I"$BUILD/mod")
SOURCES=(jfe_kinds jfe_stats jfe_performance jfe)
for source in "${SOURCES[@]}"; do
  "$FC" "${FLAGS[@]}" -c "$ROOT/src/$source.f90" -o "$BUILD/obj/$source.o"
done
OBJECTS=("$BUILD/obj/jfe_kinds.o" "$BUILD/obj/jfe_stats.o" \
  "$BUILD/obj/jfe_performance.o" "$BUILD/obj/jfe.o")
for file in "$ROOT"/test/*.f90; do
  name=$(basename "${file%.f90}")
  "$FC" "${FLAGS[@]}" "$file" "${OBJECTS[@]}" -o "$BUILD/bin/$name"
  "$BUILD/bin/$name"
done
