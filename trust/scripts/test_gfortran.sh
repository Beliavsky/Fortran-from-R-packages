#!/usr/bin/env sh
set -eu
FC=${FC:-gfortran}
FLAGS="-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0"
BUILD=build_strict
rm -rf "$BUILD"
mkdir -p "$BUILD"
for name in trust_kinds trust_types trust_linalg trust_core trust; do
    "$FC" $FLAGS -J "$BUILD" -I "$BUILD" -c "src/$name.f90" -o "$BUILD/$name.o"
done
OBJS="$BUILD/trust_kinds.o $BUILD/trust_types.o $BUILD/trust_linalg.o $BUILD/trust_core.o $BUILD/trust.o"
for source in test/*.f90; do
    exe="$BUILD/$(basename "$source" .f90)"
    "$FC" $FLAGS -J "$BUILD" -I "$BUILD" "$source" $OBJS -o "$exe"
    "$exe"
done
for source in example/*.f90; do
    exe="$BUILD/$(basename "$source" .f90)_example"
    "$FC" $FLAGS -J "$BUILD" -I "$BUILD" "$source" $OBJS -o "$exe"
    "$exe"
done
