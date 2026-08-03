#!/bin/sh
set -eu

cd "$(dirname "$0")/.."
rm -rf build-gfortran-optimized
mkdir -p build-gfortran-optimized/mod build-gfortran-optimized/obj build-gfortran-optimized/bin

flags="-std=f2018 -O3 -Wall -Wextra -Wpedantic -Werror"
sources="src/moments_kinds.f90 src/moments_status.f90 src/moments_probability.f90 src/moments_statistics.f90 src/moments_transforms.f90 src/moments_tests.f90 src/moments.f90"

for source in $sources; do
  object="build-gfortran-optimized/obj/$(basename "$source" .f90).o"
  gfortran $flags -Jbuild-gfortran-optimized/mod -Ibuild-gfortran-optimized/mod -c "$source" -o "$object"
done
ar rcs build-gfortran-optimized/libmoments.a build-gfortran-optimized/obj/*.o

for source in test/*.f90; do
  executable="build-gfortran-optimized/bin/$(basename "$source" .f90)"
  gfortran $flags -Ibuild-gfortran-optimized/mod "$source" \
    build-gfortran-optimized/libmoments.a -o "$executable"
  "$executable"
done

for source in example/*.f90 app/*.f90; do
  executable="build-gfortran-optimized/bin/$(basename "$source" .f90)"
  gfortran $flags -Ibuild-gfortran-optimized/mod "$source" \
    build-gfortran-optimized/libmoments.a -o "$executable"
  "$executable" >/dev/null
done

printf '%s\n' 'optimized GNU Fortran validation: PASS'
