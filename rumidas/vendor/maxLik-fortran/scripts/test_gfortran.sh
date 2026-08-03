#!/bin/sh
set -eu

cd "$(dirname "$0")/.."
rm -rf build-gfortran
mkdir -p build-gfortran/mod build-gfortran/obj build-gfortran/bin

flags="-std=f2018 -Wall -Wextra -Wpedantic -Werror -Wimplicit-interface -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -O0"
sources="src/maxlik_kinds.f90 src/maxlik_status.f90 src/maxlik_types.f90 src/maxlik_linalg.f90 src/maxlik_random.f90 src/maxlik_evaluation.f90 src/maxlik_solvers.f90 src/maxlik_inference.f90 src/maxlik_utilities.f90 src/maxlik_api.f90 src/maxlik.f90"

for source in $sources; do
  object="build-gfortran/obj/$(basename "$source" .f90).o"
  gfortran $flags -Jbuild-gfortran/mod -Ibuild-gfortran/mod -c "$source" -o "$object"
done
ar rcs build-gfortran/libmaxlik.a build-gfortran/obj/*.o

for source in test/*.f90; do
  executable="build-gfortran/bin/$(basename "$source" .f90)"
  gfortran $flags -Ibuild-gfortran/mod "$source" build-gfortran/libmaxlik.a -o "$executable"
  "$executable"
done

for source in example/*.f90 app/*.f90; do
  executable="build-gfortran/bin/$(basename "$source" .f90)"
  gfortran $flags -Ibuild-gfortran/mod "$source" build-gfortran/libmaxlik.a -o "$executable"
  "$executable" >/dev/null
done

printf '%s\n' 'strict GNU Fortran validation: PASS'
