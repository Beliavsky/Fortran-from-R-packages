#!/usr/bin/env sh
set -eu
fc=${FC:-gfortran}
flags='-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0'
rm -rf build-strict
mkdir build-strict
for f in src/optimflex_types.f90 src/optimflex_linalg.f90 src/optimflex_diff.f90 \
         src/optimflex_helpers.f90 src/optimflex_optimizers.f90 src/optimflex.f90; do
  $fc $flags -c -J build-strict -I build-strict "$f" -o "build-strict/$(basename "$f" .f90).o"
done
objs=$(find build-strict -name '*.o' -print | sort | tr '\n' ' ')
for f in test/*.f90; do
  exe="build-strict/$(basename "$f" .f90)"
  $fc $flags -I build-strict "$f" $objs -o "$exe"
  "$exe"
done
