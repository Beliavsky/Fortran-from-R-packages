#!/usr/bin/env sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
BUILD="$ROOT/build-strict"
rm -rf "$BUILD"
mkdir -p "$BUILD/mod" "$BUILD/obj" "$BUILD/bin"
FC=${FC:-gfortran}
FLAGS="-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all"
for src in lsei_kinds lsei_types lsei_linalg lsei_nnls lsei_solver lsei_utils lsei; do
  "$FC" $FLAGS -J "$BUILD/mod" -I "$BUILD/mod" -c "$ROOT/src/$src.f90" -o "$BUILD/obj/$src.o"
done
OBJS=$(printf ' %s' "$BUILD"/obj/*.o)
for src in "$ROOT"/test/*.f90; do
  name=$(basename "$src" .f90)
  "$FC" $FLAGS -I "$BUILD/mod" "$src" $OBJS -o "$BUILD/bin/$name"
  "$BUILD/bin/$name"
done
for src in "$ROOT"/example/*.f90; do
  name=$(basename "$src" .f90)
  "$FC" $FLAGS -I "$BUILD/mod" "$src" $OBJS -o "$BUILD/bin/$name"
  "$BUILD/bin/$name"
done
