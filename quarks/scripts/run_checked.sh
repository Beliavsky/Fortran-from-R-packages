#!/usr/bin/env sh
set -eu
cd "$(dirname "$0")/.."
. scripts/common_sources.sh
BUILD=build/checked
rm -rf "$BUILD"
mkdir -p "$BUILD/mod" "$BUILD/bin"
FLAGS="-O0 -g -std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow"
for source in $RUGARCH_SOURCES $QUARKS_SOURCES; do
  object="$BUILD/$(basename "${source%.f90}").o"
  gfortran $FLAGS -J "$BUILD/mod" -I "$BUILD/mod" -c "$source" -o "$object"
done
ar rcs "$BUILD/libquarks.a" "$BUILD"/*.o
for source in test/*.f90; do
  name=$(basename "${source%.f90}")
  gfortran $FLAGS -J "$BUILD/mod" -I "$BUILD/mod" "$source" \
    "$BUILD/libquarks.a" -o "$BUILD/bin/$name"
  "$BUILD/bin/$name"
done
for source in example/*.f90 app/*.f90; do
  name=$(basename "${source%.f90}")
  gfortran $FLAGS -J "$BUILD/mod" -I "$BUILD/mod" "$source" \
    "$BUILD/libquarks.a" -o "$BUILD/bin/$name"
  "$BUILD/bin/$name"
done
