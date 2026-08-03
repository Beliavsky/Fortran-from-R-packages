#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "$0")/.." && pwd)"
build="$root/build/gfortran-release"
rm -rf "$build"
mkdir -p "$build"
cd "$build"

flags=(-std=f2018 -O3 -Wall -Wextra -Wpedantic -Werror -J. -I.)
sources=(
  "$root/src/intraday_kinds.f90"
  "$root/src/intraday_types.f90"
  "$root/src/intraday_utils.f90"
  "$root/src/intraday_kalman.f90"
  "$root/src/intraday_fit.f90"
  "$root/src/intraday_use.f90"
  "$root/src/intraday_simulation.f90"
  "$root/src/intraday_model.f90"
)

gfortran "${flags[@]}" -c "${sources[@]}" "$root/test/test_support.f90"
objects=(intraday_kinds.o intraday_types.o intraday_utils.o intraday_kalman.o \
         intraday_fit.o intraday_use.o intraday_simulation.o intraday_model.o)

for source in "$root"/test/test_*.f90; do
  [[ "$(basename "$source")" == "test_support.f90" ]] && continue
  exe="$(basename "${source%.f90}")"
  gfortran "${flags[@]}" "${objects[@]}" test_support.o "$source" -o "$exe"
  "./$exe"
done

for source in "$root"/example/*.f90 "$root"/app/*.f90; do
  exe="$(basename "${source%.f90}")"
  gfortran "${flags[@]}" "${objects[@]}" "$source" -o "$exe"
  "./$exe"
done
