#!/usr/bin/env sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
SOURCES="jrvfinance_kinds.f90 jrvfinance_types.f90 jrvfinance_dates.f90 jrvfinance_roots.f90 jrvfinance_cashflows.f90 jrvfinance_bonds.f90 jrvfinance_options.f90 jrvfinance.f90"
run_config() {
  NAME=$1
  FLAGS=$2
  BUILD="$ROOT/build-$NAME"
  rm -rf "$BUILD"
  mkdir -p "$BUILD/mod"
  cd "$BUILD"
  for SRC in $SOURCES; do
    gfortran $FLAGS -J mod -I mod -c "$ROOT/src/$SRC"
  done
  for SRC in "$ROOT"/test/*.f90; do
    EXE=$(basename "$SRC" .f90)
    gfortran $FLAGS -J mod -I mod "$SRC" ./*.o -o "$EXE"
    "./$EXE"
  done
  for SRC in "$ROOT"/example/*.f90 "$ROOT"/app/*.f90; do
    EXE=$(basename "$SRC" .f90)
    gfortran $FLAGS -J mod -I mod "$SRC" ./*.o -o "$EXE"
    "./$EXE" >/dev/null
  done
}
run_config debug "-std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow"
run_config release "-std=f2018 -O3 -Wall -Wextra -Wpedantic -Werror"
echo "All jrvFinance-fortran validations passed."
