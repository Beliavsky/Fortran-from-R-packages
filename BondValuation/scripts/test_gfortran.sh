#!/usr/bin/env sh
# SPDX-License-Identifier: GPL-3.0-only
set -eu

mode=${1:-debug}
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="$root/build-gfortran-$mode"
rm -rf "$build"
mkdir -p "$build/obj" "$build/mod" "$build/bin"

common="-std=f2018 -Wall -Wextra -Werror -ffree-line-length-none -J$build/mod -I$build/mod"
case "$mode" in
  debug) flags="$common -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow" ;;
  release) flags="$common -O3" ;;
  *) echo "usage: $0 [debug|release]" >&2; exit 2 ;;
esac

sources="bondvaluation_kinds bondvaluation_brazil_calendar bondvaluation_dates bondvaluation_daycount bondvaluation_schedule bondvaluation_pricing bondvaluation_compat bondvaluation"
objects=""
for name in $sources; do
  gfortran $flags -c "$root/src/$name.f90" -o "$build/obj/$name.o"
  objects="$objects $build/obj/$name.o"
done

for source in "$root"/test/*.f90 "$root"/app/*.f90 "$root"/example/*.f90; do
  name=$(basename "$source" .f90)
  gfortran $flags "$source" $objects -o "$build/bin/$name"
done

for test_program in "$build"/bin/test_*; do
  "$test_program"
done

"$build/bin/bondvaluation_demo" >/dev/null
"$build/bin/day_count_comparison" >/dev/null
"$build/bin/regular_bond" >/dev/null

echo "GNU Fortran $mode build: PASS"
