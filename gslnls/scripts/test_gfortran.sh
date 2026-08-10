#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="$ROOT/build-strict"
rm -rf "$BUILD"
mkdir -p "$BUILD/mod" "$BUILD/bin"
cd "$BUILD"
FC="${FC:-gfortran}"
FLAGS="${FFLAGS:--std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all}"
SOURCES=(
  "$ROOT/src/gslnls_kinds.f90"
  "$ROOT/src/gslnls_types.f90"
  "$ROOT/src/gslnls_linalg.f90"
  "$ROOT/src/gslnls_loss.f90"
  "$ROOT/src/gslnls_core.f90"
  "$ROOT/src/gslnls_multistart.f90"
  "$ROOT/src/gslnls_stats.f90"
  "$ROOT/src/gslnls.f90"
)
for src in "${SOURCES[@]}"; do
  obj="$(basename "${src%.f90}").o"
  "$FC" $FLAGS -J"$BUILD/mod" -I"$BUILD/mod" -c "$src" -o "$obj"
done
OBJS=(gslnls_kinds.o gslnls_types.o gslnls_linalg.o gslnls_loss.o gslnls_core.o gslnls_multistart.o gslnls_stats.o gslnls.o)
for src in "$ROOT"/test/*.f90; do
  exe="$BUILD/bin/$(basename "${src%.f90}")"
  "$FC" $FLAGS -J"$BUILD/mod" -I"$BUILD/mod" "$src" "${OBJS[@]}" -o "$exe"
  "$exe"
done
for src in "$ROOT"/example/*.f90; do
  exe="$BUILD/bin/$(basename "${src%.f90}")"
  "$FC" $FLAGS -J"$BUILD/mod" -I"$BUILD/mod" "$src" "${OBJS[@]}" -o "$exe"
  "$exe"
done
