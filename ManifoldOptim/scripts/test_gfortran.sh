#!/usr/bin/env sh
set -eu
rm -rf build_gfortran
mkdir build_gfortran
FC=${FC:-gfortran}
FLAGS="-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all"
$FC $FLAGS -J build_gfortran -I build_gfortran -c \
  src/manifoldoptim_kinds.f90 \
  src/manifoldoptim_types.f90 \
  src/manifoldoptim_linalg.f90 \
  src/manifoldoptim_manifolds.f90 \
  src/manifoldoptim_solvers.f90 \
  src/manifoldoptim.f90
for t in test/*.f90; do
  exe="build_gfortran/$(basename "$t" .f90)"
  $FC $FLAGS -J build_gfortran -I build_gfortran ./*.o "$t" -o "$exe"
  "$exe"
done
rm -f ./*.o
