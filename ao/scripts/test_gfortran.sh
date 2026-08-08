#!/usr/bin/env sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
BUILD="$ROOT/build-gfortran"
rm -rf "$BUILD"
mkdir -p "$BUILD"
cd "$BUILD"
FLAGS="-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0"
gfortran $FLAGS -J. -I. -c \
  "$ROOT/src/ao_kinds.f90" \
  "$ROOT/src/ao_types.f90" \
  "$ROOT/src/ao_random.f90" \
  "$ROOT/src/ao_history.f90" \
  "$ROOT/src/ao_base_optimizer.f90" \
  "$ROOT/src/ao.f90"
for src in "$ROOT"/test/*.f90; do
  exe="$(basename "$src" .f90)"
  gfortran $FLAGS -I. ./*.o "$src" -o "$exe"
  "./$exe"
done
for src in "$ROOT"/example/*.f90; do
  exe="$(basename "$src" .f90)"
  gfortran $FLAGS -I. ./*.o "$src" -o "$exe"
  "./$exe"
done
