#!/usr/bin/env sh
set -eu
rm -rf build_strict
mkdir build_strict
FC=${FC:-gfortran}
FLAGS='-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0'
for s in src/smoof_kinds.f90 src/smoof_single.f90 src/smoof_multi.f90 \
  src/smoof_cec09.f90 src/smoof_cec2019.f90 src/smoof_ed.f90 \
  src/smoof_nk.f90 src/smoof.f90; do
  "$FC" $FLAGS -Ibuild_strict -Jbuild_strict -c "$s" \
    -o "build_strict/$(basename "${s%.f90}").o"
done
OBJS=$(printf '%s ' build_strict/*.o)
for t in test/*.f90; do
  exe="build_strict/$(basename "${t%.f90}")"
  "$FC" $FLAGS -Ibuild_strict "$t" $OBJS -o "$exe"
  "$exe"
done
for e in example/*.f90; do
  exe="build_strict/$(basename "${e%.f90}")"
  "$FC" $FLAGS -Ibuild_strict "$e" $OBJS -o "$exe"
  "$exe"
done
