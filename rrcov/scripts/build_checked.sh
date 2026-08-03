#!/usr/bin/env sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="$root/build/checked"
rm -rf "$build"
mkdir -p "$build"
cd "$build"
flags="-std=f2018 -Wall -Wextra -Werror -fcheck=all -fbacktrace -fimplicit-none -ffree-line-length-none -O0 -g"
sources="
$root/src/rrcov_kinds.f90
$root/src/rrcov_types.f90
$root/src/rrcov_random.f90
$root/src/rrcov_sort.f90
$root/src/rrcov_linalg.f90
$root/src/rrcov_stats.f90
$root/src/rrcov_robust.f90
$root/src/rrcov_pca.f90
$root/src/rrcov_da.f90
$root/src/rrcov_tests.f90
$root/src/rrcov_utils.f90
$root/src/rrcov.f90
"
gfortran $flags -c $sources
for test_source in "$root"/test/*.f90; do
  test_name=$(basename "$test_source" .f90)
  gfortran $flags "$test_source" ./*.o -o "$test_name"
  ./"$test_name"
done
gfortran $flags "$root/app/demo_rrcov.f90" ./*.o -o demo_rrcov
./demo_rrcov
gfortran $flags "$root/example/example_covariance.f90" ./*.o -o example_covariance
./example_covariance
