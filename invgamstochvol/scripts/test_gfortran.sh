#!/bin/sh
set -eu

cd "$(dirname "$0")/.."
rm -rf build-gfortran
mkdir -p build-gfortran/mod build-gfortran/obj build-gfortran/bin

flags="-std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace"
sources="src/invgamstochvol_kinds.f90 src/invgamstochvol_status.f90 src/invgamstochvol_rng.f90 src/invgamstochvol_special.f90 src/invgamstochvol_model.f90 src/invgamstochvol.f90"

for source in $sources; do
  object="build-gfortran/obj/$(basename "$source" .f90).o"
  gfortran $flags -Jbuild-gfortran/mod -Ibuild-gfortran/mod -c "$source" -o "$object"
done
ar rcs build-gfortran/libinvgamstochvol.a build-gfortran/obj/*.o

for source in test/*.f90; do
  executable="build-gfortran/bin/$(basename "$source" .f90)"
  gfortran $flags -Ibuild-gfortran/mod "$source" build-gfortran/libinvgamstochvol.a -o "$executable"
  "$executable"
done

for source in example/*.f90 app/*.f90; do
  executable="build-gfortran/bin/$(basename "$source" .f90)"
  gfortran $flags -Ibuild-gfortran/mod "$source" build-gfortran/libinvgamstochvol.a -o "$executable"
  "$executable" >/dev/null
done

printf '%s\n' 'strict GNU Fortran validation: PASS'
