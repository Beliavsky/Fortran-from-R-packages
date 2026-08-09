#!/usr/bin/env sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
BUILD="$ROOT/build-strict"
rm -rf "$BUILD"
mkdir -p "$BUILD"
FFLAGS="-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0"
gfortran $FFLAGS -J "$BUILD" -I "$BUILD" -c "$ROOT/src/lowrankqp_kinds.f90" -o "$BUILD/kinds.o"
gfortran $FFLAGS -J "$BUILD" -I "$BUILD" -c "$ROOT/src/lowrankqp_linalg.f90" -o "$BUILD/linalg.o"
gfortran $FFLAGS -J "$BUILD" -I "$BUILD" -c "$ROOT/src/lowrankqp.f90" -o "$BUILD/lowrankqp.o"
for src in "$ROOT"/test/*.f90; do
    exe="$BUILD/$(basename "${src%.f90}")"
    gfortran $FFLAGS -J "$BUILD" -I "$BUILD" "$BUILD/kinds.o" "$BUILD/linalg.o" "$BUILD/lowrankqp.o" "$src" -o "$exe"
    "$exe"
done
for src in "$ROOT"/example/*.f90; do
    exe="$BUILD/$(basename "${src%.f90}")"
    gfortran $FFLAGS -J "$BUILD" -I "$BUILD" "$BUILD/kinds.o" "$BUILD/linalg.o" "$BUILD/lowrankqp.o" "$src" -o "$exe"
    "$exe" >/dev/null
done
