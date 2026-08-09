#!/usr/bin/env sh
set -eu
FC=${FC:-gfortran}
FLAGS="-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0"
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
BUILD="$ROOT/build-strict"
rm -rf "$BUILD"
mkdir -p "$BUILD/mod" "$BUILD/obj" "$BUILD/bin"

compile () {
  src="$1"
  obj="$BUILD/obj/$(basename "${src%.f90}").o"
  "$FC" $FLAGS -I"$BUILD/mod" -J"$BUILD/mod" -c "$ROOT/src/$src" -o "$obj"
}

compile qpoases_kinds.f90
compile qpoases_types.f90
compile qpoases_linalg.f90
compile qpoases_active_set.f90
compile qpoases_solver.f90
compile qpoases.f90
compile roi_qpoases.f90

OBJS="$BUILD/obj/qpoases_kinds.o $BUILD/obj/qpoases_types.o \
$BUILD/obj/qpoases_linalg.o $BUILD/obj/qpoases_active_set.o \
$BUILD/obj/qpoases_solver.o $BUILD/obj/qpoases.o $BUILD/obj/roi_qpoases.o"

for src in "$ROOT"/test/*.f90; do
  name=$(basename "${src%.f90}")
  "$FC" $FLAGS -I"$BUILD/mod" "$src" $OBJS -o "$BUILD/bin/$name"
  "$BUILD/bin/$name"
done

for src in "$ROOT"/example/*.f90; do
  name=$(basename "${src%.f90}")
  "$FC" $FLAGS -I"$BUILD/mod" "$src" $OBJS -o "$BUILD/bin/$name"
  "$BUILD/bin/$name"
done
