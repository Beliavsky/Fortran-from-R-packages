#!/bin/sh
set -eu

cd "$(dirname "$0")/.."
rm -rf build-gfortran-optimized
mkdir -p build-gfortran-optimized/mod build-gfortran-optimized/obj build-gfortran-optimized/bin

flags="-std=f2018 -Wall -Wextra -Wpedantic -Werror -Wimplicit-interface -O3"
maxlik_sources="vendor/maxLik-fortran/src/maxlik_kinds.f90 vendor/maxLik-fortran/src/maxlik_status.f90 vendor/maxLik-fortran/src/maxlik_types.f90 vendor/maxLik-fortran/src/maxlik_linalg.f90 vendor/maxLik-fortran/src/maxlik_random.f90 vendor/maxLik-fortran/src/maxlik_evaluation.f90 vendor/maxLik-fortran/src/maxlik_solvers.f90 vendor/maxLik-fortran/src/maxlik_inference.f90 vendor/maxLik-fortran/src/maxlik_utilities.f90 vendor/maxLik-fortran/src/maxlik_api.f90 vendor/maxLik-fortran/src/maxlik.f90"
rumidas_sources="src/rumidas_kinds.f90 src/rumidas_status.f90 src/rumidas_types.f90 src/rumidas_weights.f90 src/rumidas_statistics.f90 src/rumidas_garch_midas.f90 src/rumidas_mem.f90 src/rumidas_fit.f90 src/rumidas_forecast.f90 src/rumidas.f90"

for source in $maxlik_sources $rumidas_sources; do
  object="build-gfortran-optimized/obj/$(basename "$source" .f90).o"
  gfortran $flags -Jbuild-gfortran-optimized/mod -Ibuild-gfortran-optimized/mod -c "$source" -o "$object"
done
ar rcs build-gfortran-optimized/librumidas.a build-gfortran-optimized/obj/*.o

for source in test/*.f90; do
  executable="build-gfortran-optimized/bin/$(basename "$source" .f90)"
  gfortran $flags -Ibuild-gfortran-optimized/mod "$source" build-gfortran-optimized/librumidas.a -o "$executable"
  "$executable"
done

for source in example/*.f90 app/*.f90; do
  executable="build-gfortran-optimized/bin/$(basename "$source" .f90)"
  gfortran $flags -Ibuild-gfortran-optimized/mod "$source" build-gfortran-optimized/librumidas.a -o "$executable"
  "$executable" >/dev/null
done

printf '%s\n' 'optimized GNU Fortran validation: PASS'
