#!/usr/bin/env sh
set -eu

FC=${FC:-gfortran}
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
BUILD="$ROOT/build_gfortran"
MOD="$BUILD/mod"
OBJ="$BUILD/obj"
BIN="$BUILD/bin"

rm -rf "$BUILD"
mkdir -p "$MOD" "$OBJ" "$BIN"

FLAGS="-std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow -O0 -g"
SOURCES="garchito_kinds garchito_callbacks garchito_types garchito_utils garchito_optimizer garchito_models garchito"
OBJECTS=""
for unit in $SOURCES; do
  "$FC" $FLAGS -J"$MOD" -I"$MOD" -c "$ROOT/src/$unit.f90" -o "$OBJ/$unit.o"
  OBJECTS="$OBJECTS $OBJ/$unit.o"
done

for source in "$ROOT"/test/*.f90; do
  name=$(basename "$source" .f90)
  "$FC" $FLAGS -J"$MOD" -I"$MOD" $OBJECTS "$source" -o "$BIN/$name"
  "$BIN/$name"
done

for source in "$ROOT"/example/*.f90 "$ROOT"/app/*.f90; do
  name=$(basename "$source" .f90)
  "$FC" $FLAGS -J"$MOD" -I"$MOD" $OBJECTS "$source" -o "$BIN/$name"
  "$BIN/$name"
done
