#!/usr/bin/env sh
set -eu
rm -rf build_gfortran
mkdir build_gfortran
cd build_gfortran
flags="-std=f2018 -Wall -Wextra -Werror -fcheck=all -fbacktrace"
sources="../src/fitheavytail_kinds.f90 ../src/fitheavytail_status.f90 ../src/fitheavytail_types.f90 ../src/fitheavytail_linalg.f90 ../src/fitheavytail_special.f90 ../src/fitheavytail_rng.f90 ../src/fitheavytail_tail.f90 ../src/fitheavytail_elliptical.f90 ../src/fitheavytail_mvt.f90 ../src/fitheavytail_mvst.f90 ../src/fitheavytail.f90 ../test/test_support.f90"
gfortran $flags -J. -I. -c $sources
for source in ../test/test_tail_estimators.f90 ../test/test_elliptical.f90 ../test/test_mvt.f90 ../test/test_mvst.f90; do
  name=$(basename "$source" .f90)
  gfortran $flags *.o "$source" -o "$name"
  ./$name
done
