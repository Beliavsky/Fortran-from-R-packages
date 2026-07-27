#!/usr/bin/env sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="$root/build-validation"
rm -rf "$build"
mkdir -p "$build"
cd "$build"
flags="-std=f2018 -Wall -Wextra -Wpedantic -Wconversion-extra -Wimplicit-interface -Werror -fcheck=all -fbacktrace -O0"
sources="../src/rnd_kinds.f90 ../src/rnd_types.f90 ../src/rnd_special.f90 ../src/rnd_linalg.f90 ../src/rnd_densities.f90 ../src/rnd_pricing.f90 ../src/rnd_optimize.f90 ../src/rnd_objectives.f90 ../src/rnd_fitting.f90 ../src/rnd.f90"
gfortran $flags -J . -I . -c $sources
for source in ../test/*.f90; do
  name=$(basename "$source" .f90)
  gfortran $flags -J . -I . "$source" ./*.o -o "$name"
  "./$name"
done
for source in ../app/*.f90 ../example/*.f90; do
  name=$(basename "$source" .f90)
  gfortran $flags -J . -I . "$source" ./*.o -o "$name"
  "./$name"
done
