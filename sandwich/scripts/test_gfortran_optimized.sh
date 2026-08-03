#!/bin/sh
set -eu

cd "$(dirname "$0")/.."
rm -rf build-gfortran-optimized
mkdir -p build-gfortran-optimized/mod build-gfortran-optimized/obj build-gfortran-optimized/bin

flags="-std=f2018 -O3 -Wall -Wextra -Wpedantic -Werror"
sources="src/sandwich_kinds.f90 src/sandwich_status.f90 src/sandwich_utils.f90 src/sandwich_linalg.f90 src/sandwich_core.f90 src/sandwich_regression.f90 src/sandwich_kernels.f90 src/sandwich_auxiliary.f90 src/sandwich_hc.f90 src/sandwich_hac.f90 src/sandwich_cluster.f90 src/sandwich_panel.f90 src/sandwich_bootstrap.f90 src/sandwich.f90"

for source in $sources; do
  object="build-gfortran-optimized/obj/$(basename "$source" .f90).o"
  gfortran $flags -Jbuild-gfortran-optimized/mod -Ibuild-gfortran-optimized/mod -c "$source" -o "$object"
done
ar rcs build-gfortran-optimized/libsandwich.a build-gfortran-optimized/obj/*.o

for source in test/*.f90; do
  executable="build-gfortran-optimized/bin/$(basename "$source" .f90)"
  gfortran $flags -Ibuild-gfortran-optimized/mod "$source" \
    build-gfortran-optimized/libsandwich.a -o "$executable"
  "$executable"
done

for source in example/*.f90 app/*.f90; do
  executable="build-gfortran-optimized/bin/$(basename "$source" .f90)"
  gfortran $flags -Ibuild-gfortran-optimized/mod "$source" \
    build-gfortran-optimized/libsandwich.a -o "$executable"
  "$executable" >/dev/null
done

printf '%s\n' 'optimized GNU Fortran validation: PASS'
