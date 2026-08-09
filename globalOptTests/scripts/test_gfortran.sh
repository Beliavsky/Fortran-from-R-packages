#!/usr/bin/env sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
BUILD="${TMPDIR:-/tmp}/globalOptTests-fortran-build"
rm -rf "$BUILD"
mkdir -p "$BUILD/mod" "$BUILD/bin"
FLAGS="-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0"
gfortran $FLAGS -c "$ROOT/src/global_opt_tests.f90" -J "$BUILD/mod" -o "$BUILD/global_opt_tests.o"
for src in "$ROOT"/test/*.f90; do
    exe="$BUILD/bin/$(basename "${src%.f90}")"
    gfortran $FLAGS -I "$BUILD/mod" "$src" "$BUILD/global_opt_tests.o" -o "$exe"
    "$exe"
done
for src in "$ROOT"/example/*.f90; do
    exe="$BUILD/bin/$(basename "${src%.f90}")"
    gfortran $FLAGS -I "$BUILD/mod" "$src" "$BUILD/global_opt_tests.o" -o "$exe"
    "$exe"
done
