#!/usr/bin/env sh
set -eu
FC=${FC:-gfortran}
FLAGS="-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all"
BUILD=.strict-build
rm -rf "$BUILD"
mkdir -p "$BUILD"
for f in src/oor_kinds.f90 src/oor_interfaces.f90 src/oor_random.f90 src/oor_test_functions.f90 src/oor_poo.f90 src/oor_stosoo.f90 src/oor.f90; do
  "$FC" $FLAGS -c -J "$BUILD" -I "$BUILD" "$f" -o "$BUILD/$(basename "${f%.f90}").o"
done
for t in test/*.f90; do
  exe="$BUILD/$(basename "${t%.f90}")"
  "$FC" $FLAGS -J "$BUILD" -I "$BUILD" "$t" "$BUILD"/*.o -o "$exe"
  "$exe"
done
