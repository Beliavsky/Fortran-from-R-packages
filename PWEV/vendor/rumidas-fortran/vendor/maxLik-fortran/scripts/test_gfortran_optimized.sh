#!/bin/sh
set -eu

cd "$(dirname "$0")/.."
rm -rf build-gfortran-optimized
mkdir -p build-gfortran-optimized/mod build-gfortran-optimized/obj build-gfortran-optimized/bin

flags="-std=f2018 -Wall -Wextra -Wpedantic -Werror -Wimplicit-interface -O3 -march=native"
sources="src/maxlik_kinds.f90 src/maxlik_status.f90 src/maxlik_types.f90 src/maxlik_linalg.f90 src/maxlik_random.f90 src/maxlik_evaluation.f90 src/maxlik_solvers.f90 src/maxlik_inference.f90 src/maxlik_utilities.f90 src/maxlik_api.f90 src/maxlik.f90"

for source in $sources; do
  object="build-gfortran-optimized/obj/$(basename "$source" .f90).o"
  gfortran $flags -Jbuild-gfortran-optimized/mod -Ibuild-gfortran-optimized/mod -c "$source" -o "$object"
done
ar rcs build-gfortran-optimized/libmaxlik.a build-gfortran-optimized/obj/*.o

for source in test/*.f90; do
  executable="build-gfortran-optimized/bin/$(basename "$source" .f90)"
  gfortran $flags -Ibuild-gfortran-optimized/mod "$source" build-gfortran-optimized/libmaxlik.a -o "$executable"
  "$executable"
done

for source in example/*.f90 app/*.f90; do
  executable="build-gfortran-optimized/bin/$(basename "$source" .f90)"
  gfortran $flags -Ibuild-gfortran-optimized/mod "$source" build-gfortran-optimized/libmaxlik.a -o "$executable"
  "$executable" >/dev/null
done

printf '%s\n' 'optimized GNU Fortran validation: PASS'
