#!/bin/sh
set -eu

cd "$(dirname "$0")/.."
rm -rf build-gfortran-opt
mkdir -p build-gfortran-opt/mod build-gfortran-opt/obj build-gfortran-opt/bin

flags="-O3 -std=f2018 -Wall -Wextra -Wpedantic -Werror"
sources="src/invgamstochvol_kinds.f90 src/invgamstochvol_status.f90 src/invgamstochvol_rng.f90 src/invgamstochvol_special.f90 src/invgamstochvol_model.f90 src/invgamstochvol.f90"

for source in $sources; do
  object="build-gfortran-opt/obj/$(basename "$source" .f90).o"
  gfortran $flags -Jbuild-gfortran-opt/mod -Ibuild-gfortran-opt/mod -c "$source" -o "$object"
done
ar rcs build-gfortran-opt/libinvgamstochvol.a build-gfortran-opt/obj/*.o

for source in test/*.f90; do
  executable="build-gfortran-opt/bin/$(basename "$source" .f90)"
  gfortran $flags -Ibuild-gfortran-opt/mod "$source" build-gfortran-opt/libinvgamstochvol.a -o "$executable"
  "$executable"
done

for source in example/*.f90 app/*.f90; do
  executable="build-gfortran-opt/bin/$(basename "$source" .f90)"
  gfortran $flags -Ibuild-gfortran-opt/mod "$source" build-gfortran-opt/libinvgamstochvol.a -o "$executable"
  "$executable" >/dev/null
done

printf '%s\n' 'optimized GNU Fortran validation: PASS'
