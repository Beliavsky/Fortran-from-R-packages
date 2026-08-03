#!/usr/bin/env sh
set -eu

FC=${FC:-gfortran}
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
BUILD="$ROOT/build-gfortran"
rm -rf "$BUILD"
mkdir -p "$BUILD/mod" "$BUILD/obj" "$BUILD/bin"

FLAGS=${FFLAGS:-"-std=f2018 -Wall -Wextra -Werror -pedantic -fcheck=all -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -O0 -g"}
SOURCES="fints_kinds fints_status fints_types fints_linalg fints_special fints_summary fints_time_series fints_finance fints_dates fints_apca fints_arma fints_arima fints"
OBJECTS=""
for unit in $SOURCES; do
  "$FC" $FLAGS -J "$BUILD/mod" -I "$BUILD/mod" -c "$ROOT/src/$unit.f90" -o "$BUILD/obj/$unit.o"
  OBJECTS="$OBJECTS $BUILD/obj/$unit.o"
done

for source in "$ROOT"/test/*.f90; do
  name=$(basename "$source" .f90)
  "$FC" $FLAGS -J "$BUILD/mod" -I "$BUILD/mod" "$source" $OBJECTS -o "$BUILD/bin/$name"
  "$BUILD/bin/$name"
done

for source in "$ROOT"/app/*.f90 "$ROOT"/example/*.f90; do
  name=$(basename "$source" .f90)
  "$FC" $FLAGS -J "$BUILD/mod" -I "$BUILD/mod" "$source" $OBJECTS -o "$BUILD/bin/$name"
  "$BUILD/bin/$name" >/dev/null
done

echo "All FinTS-fortran tests and runnable targets passed."
