#!/usr/bin/env sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="$root/build/gfortran-check"
rm -rf "$build"
mkdir -p "$build"

sources="$root/src/splines_kinds.f90 \
$root/src/splines_linalg.f90 \
$root/src/splines_core.f90 \
$root/src/splines_basis.f90 \
$root/src/splines.f90"

# shellcheck disable=SC2086
gfortran -std=f2018 -O0 -g -Wall -Wextra -Wpedantic -Werror \
  -fcheck=all -fbacktrace -J "$build" -I "$build" \
  -o "$build/test_splines" $sources "$root/test/test_splines.f90"
"$build/test_splines"

# shellcheck disable=SC2086
gfortran -std=f2018 -O3 -Wall -Wextra -Wpedantic -Werror \
  -J "$build" -I "$build" \
  -o "$build/test_splines_optimized" $sources "$root/test/test_splines.f90"
"$build/test_splines_optimized"

# shellcheck disable=SC2086
gfortran -std=f2018 -O2 -Wall -Wextra -Wpedantic -Werror \
  -J "$build" -I "$build" \
  -o "$build/basic_splines" $sources "$root/example/basic_splines.f90"
"$build/basic_splines"
