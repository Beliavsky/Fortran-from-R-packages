#!/usr/bin/env sh
set -eu
FC=${FC:-gfortran}
FLAGS="-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0"
rm -rf build-strict
mkdir build-strict
$FC $FLAGS -J build-strict -I build-strict -c src/subplex.f90 -o build-strict/subplex.o
for src in test/*.f90; do
    exe="build-strict/$(basename "$src" .f90)"
    $FC $FLAGS -J build-strict -I build-strict "$src" build-strict/subplex.o -o "$exe"
    "$exe"
done
for src in example/*.f90; do
    exe="build-strict/$(basename "$src" .f90)"
    $FC $FLAGS -J build-strict -I build-strict "$src" build-strict/subplex.o -o "$exe"
    "$exe"
done
