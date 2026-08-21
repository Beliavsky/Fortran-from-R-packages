#!/usr/bin/env sh
set -eu
mkdir -p build
SRC="src/boot_kinds.f90 src/boot_special.f90 src/boot_linalg.f90 src/boot_statistics.f90 \
src/boot_resampling.f90 src/boot_ci.f90 src/boot_influence.f90 src/boot_tilt.f90 \
src/boot_importance.f90 src/boot_core.f90 src/boot_timeseries.f90 src/boot_profiles.f90 \
src/boot_simplex.f90 src/boot_saddle.f90 src/boot_envelope.f90 src/boot_smoothing.f90 \
src/boot_censored.f90 src/boot_glm_diag.f90 src/boot_nested.f90 src/boot.f90"
FLAGS="-std=f2018 -O2 -Wall -Wextra -Werror -fcheck=all -J build -I build"
rm -f ./*.o test/*.exe example/*.exe
gfortran $FLAGS -c $SRC
for t in test/*.f90; do
    e="${t%.f90}.exe"
    gfortran $FLAGS "$t" ./*.o -o "$e"
    "$e"
done
for x in example/*.f90; do
    e="${x%.f90}.exe"
    gfortran $FLAGS "$x" ./*.o -o "$e"
    "$e"
done
