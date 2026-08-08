#!/usr/bin/env sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
BUILD="$ROOT/build_strict"
rm -rf "$BUILD"
mkdir -p "$BUILD"
FFLAGS="-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all"
SRC="calibrar_kinds calibrar_interfaces calibrar_utils calibrar_random calibrar_gradient calibrar_fitness calibrar_splines calibrar_objective calibrar_stopping calibrar_test_functions calibrar_optimization calibrar"
OBJS=""
for name in $SRC; do
  gfortran $FFLAGS -c -J "$BUILD" -I "$BUILD" "$ROOT/src/$name.f90" -o "$BUILD/$name.o"
  OBJS="$OBJS $BUILD/$name.o"
done
for file in "$ROOT"/test/*.f90; do
  name=$(basename "$file" .f90)
  gfortran $FFLAGS -I "$BUILD" "$file" $OBJS -o "$BUILD/$name"
  "$BUILD/$name"
done
for file in "$ROOT"/example/*.f90; do
  name=$(basename "$file" .f90)
  gfortran $FFLAGS -I "$BUILD" "$file" $OBJS -o "$BUILD/$name"
  "$BUILD/$name"
done
