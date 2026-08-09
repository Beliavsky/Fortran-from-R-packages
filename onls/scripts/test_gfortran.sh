#!/usr/bin/env sh
set -eu
FC=${FC:-gfortran}
FFLAGS='-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0'
BUILD=${TMPDIR:-/tmp}/onls-fortran-strict
rm -rf "$BUILD"
mkdir -p "$BUILD/mod" "$BUILD/obj" "$BUILD/bin"
for f in src/onls_kinds.f90 src/onls_linalg.f90 src/onls_minimize.f90 src/onls_lm.f90 src/onls_core.f90 src/onls.f90; do
  b=$(basename "$f" .f90)
  $FC $FFLAGS -J "$BUILD/mod" -I "$BUILD/mod" -c "$f" -o "$BUILD/obj/$b.o"
done
OBJS="$BUILD/obj/onls_kinds.o $BUILD/obj/onls_linalg.o $BUILD/obj/onls_minimize.o $BUILD/obj/onls_lm.o $BUILD/obj/onls_core.o $BUILD/obj/onls.o"
for f in test/*.f90; do
  b=$(basename "$f" .f90)
  $FC $FFLAGS -I "$BUILD/mod" "$f" $OBJS -o "$BUILD/bin/$b"
  "$BUILD/bin/$b"
done
for f in example/*.f90; do
  b=$(basename "$f" .f90)
  $FC $FFLAGS -I "$BUILD/mod" "$f" $OBJS -o "$BUILD/bin/$b"
  "$BUILD/bin/$b"
done
