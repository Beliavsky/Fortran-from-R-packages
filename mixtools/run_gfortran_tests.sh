#!/usr/bin/env sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
build="$root/build-gfortran"
rm -rf "$build"
mkdir -p "$build"
cd "$build"
flags="-std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -O0"
sources="
$root/src/mixtools_kinds.f90
$root/src/mixtools_status.f90
$root/src/mixtools_types.f90
$root/src/mixtools_linalg.f90
$root/src/mixtools_rng.f90
$root/src/mixtools_distributions.f90
$root/src/mixtools_utilities.f90
$root/src/mixtools_parametric.f90
$root/src/mixtools_regression.f90
$root/src/mixtools_semiparametric.f90
$root/src/mixtools_reliability.f90
$root/src/mixtools_support.f90
$root/src/mixtools_diagnostics.f90
$root/src/mixtools_compat.f90
$root/src/mixtools.f90"
# shellcheck disable=SC2086
gfortran $flags -c $sources
for test_source in "$root"/test/*.f90; do
  test_name=$(basename "$test_source" .f90)
  # shellcheck disable=SC2086
  gfortran $flags "$test_source" ./*.o -o "$test_name"
  "./$test_name"
done
